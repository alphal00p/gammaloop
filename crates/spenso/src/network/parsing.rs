use symbolica::atom::{AtomCore, FunctionBuilder, Symbol};
use symbolica::domains::rational::Rational;
use symbolica::printer::PrintState;
use symbolica::{symbol, tag};

use super::*;

use library::Library;

use crate::network::library::DummyLibrary;
use crate::network::library::function_lib::Wrap;
use crate::network::library::panicing::ErroringLibrary;
use crate::shadowing::symbolica_utils::{AtomCoreExt, SpensoPrintSettings};
use crate::structure::abstract_index::AIND_SYMBOLS;
use crate::structure::representation::Representation;
// use crate::shadowing::Concretize;
use crate::structure::slot::{DummyAind, ParseableAind, Slot, SlotError};
use crate::structure::{
    NamedStructure, OrderedStructure, PermutedStructure, StructureContract, StructureError,
    TensorShell,
};
use crate::tensors::parametric::ParamTensor;
use crate::tensors::symbolic::SymbolicTensor;

use std::fmt::Display;

use store::TensorScalarStore;
// use log::trace;

use symbolica::atom::{AddView, Atom, AtomView, MulView, PowView, representation::FunView};

use crate::structure::{HasStructure, TensorStructure};

use crate::{shadowing::Concretize, structure::HasName, structure::representation::LibraryRep};

pub struct SpensoTags {
    pub tag: String,
    pub upper: String,
    pub lower: String,
    pub bracket: Symbol,
    pub pure_scalar: Symbol,
    pub dot: Symbol,
}

pub static SPENSO_TAG: std::sync::LazyLock<SpensoTags> = std::sync::LazyLock::new(|| SpensoTags {
    tag: tag!("broadcast"),
    upper: tag!("upper"),
    lower: tag!("lower"),
    bracket: symbol!("bracket"),
    pure_scalar: symbol!("pure_scalar"),
    dot: symbol!("dot";Symmetric,Linear; print = |a,opt|{
    match opt.custom_print_mode {
        Some(("spenso",i))=>{
            let SpensoPrintSettings{
                parens,
                commas,..
            } = SpensoPrintSettings::from(i);


            let AtomView::Fun(f) = a else {
                return None;
            };

            if f.get_nargs() != 3 {
                return None;
            }
            let mut args = f.iter();

            let a = args.next().unwrap();
            let b = args.next().unwrap();
            let c = args.next().unwrap();

            fn is_rep(view:AtomView<'_>)->bool{
                match view {
                    AtomView::Fun(f) if f.get_symbol().has_tag(&SPENSO_TAG.upper) => true,
                    AtomView::Var(s) if s.get_symbol().has_tag(&SPENSO_TAG.upper) => true,
                    _=>false
                }
            }

            let (a,b,c) = if is_rep(a) && !is_rep(b) && !is_rep(c) {
                (a,b,c)
            } else if is_rep(b) && !is_rep(a) && !is_rep(c) {
                (b,c,a)
            } else if is_rep(c) && !is_rep(a) && !is_rep(b) {
                (c,a,b)
            } else { return None};

            let mut s = String::new();
            if parens {
                s.push('(');
            }
            b.format(&mut s, opt,PrintState::new()).unwrap();
            s.push('.');
            c.format(&mut s, opt,PrintState::new()).unwrap();
            if parens {
                s.push(')');
            }
            Some(s)

        },
        _=>None
    }



    }),
});

pub type ShadowedStructure<Aind> = NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind>;

impl<'a, Aind: ParseableAind + AbsInd> TryFrom<AtomView<'a>>
    for PermutedStructure<ShadowedStructure<Aind>>
{
    type Error = StructureError;
    fn try_from(value: AtomView<'a>) -> Result<Self, Self::Error> {
        if let AtomView::Fun(f) = value {
            f.try_into()
        } else {
            Err(StructureError::ParsingError(value.to_plain_string()))
        }
    }
}

impl<'a, Aind: ParseableAind + AbsInd> TryFrom<AtomView<'a>>
    for PermutedStructure<SymbolicTensor<Aind>>
{
    type Error = StructureError;
    fn try_from(value: AtomView<'a>) -> Result<Self, Self::Error> {
        let structure = OrderedStructure::from_atomcore_unchecked(value)?;
        Ok(PermutedStructure::identity(SymbolicTensor {
            structure,
            expression: value.to_owned(),
        }))
    }
}

impl<Aind: ParseableAind + AbsInd> OrderedStructure<LibraryRep, Aind> {
    pub fn from_atomcore_unchecked<A: AtomCore>(value: A) -> Result<Self, StructureError> {
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
                    let s = Self::from_atomcore_unchecked(i)?;
                    structure = structure.merge(&s)?.0
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
                    // println!("Fun{}", f.as_view());
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
                                is_structure = Some(e.into());
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
        value.as_view().try_into()
    }
}

impl<'a, Aind: ParseableAind + AbsInd> TryFrom<&'a Atom>
    for PermutedStructure<ShadowedStructure<Aind>>
{
    type Error = StructureError;
    fn try_from(value: &'a Atom) -> Result<Self, Self::Error> {
        value.as_view().try_into()
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
    pub parse_inner_products: bool,
}

impl Default for ParseSettings {
    fn default() -> Self {
        Self {
            precontract_scalars: true,
            take_first_term_from_sum: false,
            depth_limit: None,
            parse_inner_products: true,
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
        S: TensorStructure + Clone,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        PermutedStructure<S>: for<'r> TryFrom<AtomView<'r>>,
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
        S: TensorStructure + Clone,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        PermutedStructure<S>: for<'r> TryFrom<AtomView<'r>>,
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
    fn as_leaf<S>(value: AtomView<'a>) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        PermutedStructure<S>: for<'r> TryFrom<AtomView<'r>>,
    {
        let s: Result<PermutedStructure<S>, _> = value.try_into();

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
        S: TensorStructure + Clone,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        PermutedStructure<S>: for<'r> TryFrom<AtomView<'r>>,
    {
        if let Some(a) = settings.depth_limit
            && a < state.depth
        {
            return Self::as_leaf::<S>(value.as_view());
        }

        state.depth += 1;
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

            Ok(first.n_add(rest?))
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
        S: TensorStructure + Clone,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        PermutedStructure<S>: for<'r> TryFrom<AtomView<'r>>,
        // <PermutedStructure<S>>::Error: Debug,
    {
        let symbol = value.get_symbol();

        if symbol == SPENSO_TAG.bracket {
            let mut n_muls = value
                .iter()
                .map(|a| Self::try_from_view_impl(a, state, library, settings))
                .collect::<Result<Vec<_>, _>>()?;

            Ok(n_muls.pop().unwrap().n_mul(n_muls))
        } else if symbol == SPENSO_TAG.pure_scalar {
            if value.get_nargs() != 1 {
                return Err(TensorNetworkError::TooManyArgsFunction(
                    value.as_view().to_plain_string(),
                ));
            }

            Ok(Self::from_scalar(value.iter().next().unwrap().try_into()?))
        } else if symbol == SPENSO_TAG.dot && settings.parse_inner_products {
            if value.get_nargs() != 3 {
                return Err(TensorNetworkError::TooManyArgsFunction(
                    value.as_view().to_plain_string(),
                ));
            }

            let mut args = value.iter();
            let a = args.next().unwrap();
            let b = args.next().unwrap();
            let c = args.next().unwrap();

            let (rep, l, r) = if let Ok(rep) = Representation::<LibraryRep>::try_from(a)
                && Representation::<LibraryRep>::try_from(b).is_err()
                && Representation::<LibraryRep>::try_from(c).is_err()
            {
                (rep, b, c)
            } else if let Ok(rep) = Representation::<LibraryRep>::try_from(b)
                && Representation::<LibraryRep>::try_from(a).is_err()
                && Representation::<LibraryRep>::try_from(c).is_err()
            {
                (rep, a, c)
            } else if let Ok(rep) = Representation::<LibraryRep>::try_from(c)
                && Representation::<LibraryRep>::try_from(b).is_err()
                && Representation::<LibraryRep>::try_from(a).is_err()
            {
                (rep, a, b)
            } else {
                return Err(TensorNetworkError::InvalidDotFunction(
                    value.as_view().to_plain_string(),
                ));
            };

            let dummy_index: Slot<LibraryRep> = rep.slot(1);

            fn index<D: Display>(
                ind: Slot<LibraryRep>,
                view: AtomView,
            ) -> Result<Atom, TensorNetworkError<D, Symbol>> {
                match view {
                    AtomView::Fun(f) => {
                        let mut fb = FunctionBuilder::new(f.get_symbol());

                        for a in f.iter() {
                            fb = fb.add_arg(a);
                        }
                        Ok(fb.add_arg(ind.to_atom()).finish())
                    }
                    AtomView::Var(s) => Ok(FunctionBuilder::new(s.get_symbol())
                        .add_arg(ind.to_atom())
                        .finish()),
                    a => Err(TensorNetworkError::InvalidDotFunction(a.to_plain_string())),
                }
            }

            let indexed_l = index(dummy_index, l)?;
            // let Ok(indexed_l_structure): Result<PermutedStructure<S>, _> =
            //     indexed_l.as_view().try_into()
            // else {
            //     return Err(TensorNetworkError::InvalidDotFunction(
            //         indexed_l.as_view().to_plain_string(),
            //     ));
            // };
            let indexed_r = index(dummy_index, r)?;

            Self::try_from_view_impl((indexed_l * indexed_r).as_view(), state, library, settings)
            // let Ok(indexed_r_structure): Result<PermutedStructure<S>, _> =
            //     indexed_r.as_view().try_into()
            // else {
            //     return Err(TensorNetworkError::InvalidDotFunction(
            //         indexed_l.as_view().to_plain_string(),
            //     ));
            // };

            // let l = match library.key_for_structure(&indexed_l_structure) {
            //     Ok(key) => {
            //         // println!("Adding lib");
            //         // let t = library.get(&key).unwrap();
            //         Self::library_tensor(
            //             &indexed_l_structure.structure,
            //             PermutedStructure {
            //                 structure: key,
            //                 rep_permutation: indexed_l_structure.rep_permutation,
            //                 index_permutation: indexed_l_structure.index_permutation,
            //             },
            //         )
            //     }
            //     Err(_) => Self::from_tensor(
            //         indexed_l_structure
            //             .structure
            //             .to_shell()
            //             .concretize(Some(indexed_l_structure.index_permutation.inverse())),
            //     ),
            // };

            // let r = match library.key_for_structure(&indexed_r_structure) {
            //     Ok(key) => {
            //         // println!("Adding lib");
            //         // let t = library.get(&key).unwrap();
            //         Self::library_tensor(
            //             &indexed_r_structure.structure,
            //             PermutedStructure {
            //                 structure: key,
            //                 rep_permutation: indexed_r_structure.rep_permutation,
            //                 index_permutation: indexed_r_structure.index_permutation,
            //             },
            //         )
            //     }
            //     Err(_) => Self::from_tensor(
            //         indexed_r_structure
            //             .structure
            //             .to_shell()
            //             .concretize(Some(indexed_r_structure.index_permutation.inverse())),
            //     ),
            // };

            // let res
            // Ok(l * r)
            // Ok(Self::from_scalar(a.try_into()?))
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
            let s: Result<PermutedStructure<S>, _> = value.as_view().try_into();

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
    fn try_from_pow<S, Lib: Library<S, Key = K>>(
        value: PowView<'a>,
        mut state: ParseState,
        library: &Lib,
        settings: &ParseSettings,
    ) -> std::result::Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        PermutedStructure<S>: for<'r> TryFrom<AtomView<'r>>,
    {
        if let Some(a) = settings.depth_limit
            && a < state.depth
        {
            return Self::as_leaf::<S>(value.as_view());
        }
        state.depth += 1;
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
        S: TensorStructure + Clone,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        PermutedStructure<S>: for<'r> TryFrom<AtomView<'r>>,
    {
        if let Some(a) = settings.depth_limit
            && a < state.depth
        {
            return Self::as_leaf::<S>(value.as_view());
        }
        state.depth += 1;
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

pub type SymbolicNet<Aind> =
    Network<NetworkStore<SymbolicTensor<Aind>, Atom>, DummyKey, Symbol, Aind>;

// pub type HepTensor<Aind> = MixedTensor<f64, ShadowedStructure<Aind>>;

pub type ParamNet<Aind> =
    Network<NetworkStore<ParamTensor<ShadowedStructure<Aind>>, Atom>, DummyKey, Symbol, Aind>;
// pub type ParamNet<Aind> =
//     Network<NetworkStore<ParamTensor<SymbolicTensor<Aind>, Atom>, DummyKey, Symbol, Aind>;

impl<Aind: AbsInd + DummyAind + ParseableAind + 'static> SymbolicNet<Aind> {
    pub fn snapshot_dot(&self) -> String {
        self.dot_display_impl(
            |a| a.to_bare_ordered_string(),
            |_| None,
            |a| a.expression.to_bare_ordered_string(),
            |a| a.to_string(),
        )
    }

    pub fn simple_execute(mut self) -> Atom {
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

        self.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &Wrap {})
            .unwrap();

        match self.result().unwrap() {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(a) => match a {
                TensorOrScalarOrKey::Key { .. } => panic!("aaa"),
                TensorOrScalarOrKey::Scalar(s) => s.clone(),
                TensorOrScalarOrKey::Tensor { tensor, .. } => tensor.expression.clone(),
            },
        }
    }
}

impl<Aind: AbsInd + DummyAind + ParseableAind + 'static> ParamNet<Aind> {
    pub fn simple_execute(&mut self) {
        let lib = DummyLibrary::<_>::new();

        self.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &ErroringLibrary::new())
            .unwrap();
    }
}
pub trait NetworkParse {
    #[allow(clippy::result_large_err)]
    fn parse_to_symbolic_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>>;

    #[allow(clippy::result_large_err)]
    fn parse_to_atom_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>>;
}

impl NetworkParse for Atom {
    fn parse_to_symbolic_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        self.as_view().parse_to_symbolic_net::<Aind>(settings)
    }

    fn parse_to_atom_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        self.as_view().parse_to_atom_net::<Aind>(settings)
    }
}

impl NetworkParse for AtomView<'_> {
    fn parse_to_symbolic_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<SymbolicNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

        SymbolicNet::<Aind>::try_from_view::<SymbolicTensor<Aind>, _>(*self, &lib, settings)
    }

    fn parse_to_atom_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        let lib = DummyLibrary::<ParamTensor<ShadowedStructure<Aind>>>::new();

        ParamNet::<Aind>::try_from_view::<ShadowedStructure<Aind>, _>(*self, &lib, settings)
    }
}

#[cfg(test)]
pub mod test {
    macro_rules! symbol_set {
        // Identifier symbols with explicit struct and static names
        ($struct_name:ident, $static_name:ident; $($char:ident)*) => {
            #[allow(non_snake_case)]
            pub struct $struct_name {
                $(pub $char: Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($char: symbol!(stringify!($char)),)*
            });
        };

        // String literals as individual statics
        (statics; $($field:ident : $string:literal),* $(,)?) => {
            $(
                #[allow(non_upper_case_globals)]
                pub static $field: ::std::sync::LazyLock<Symbol> = ::std::sync::LazyLock::new(|| symbol!($string));
            )*
        };

        // String literals grouped in a struct
        ($struct_name:ident, $static_name:ident; $($field:ident : $string:literal),* $(,)?) => {
            #[allow(non_snake_case)]
            pub struct $struct_name {
                $(pub $field: Symbol,)*
            }

            pub static $static_name: ::std::sync::LazyLock<$struct_name> = ::std::sync::LazyLock::new(|| $struct_name {
                $($field: symbol!($string),)*
            });
        };
    }

    // Generate TestSymbols with all alphabet characters and some multi-character symbols
    symbol_set!(TestSymbols, TS;
        a b c d e f h i j k l m n o p q r s t u v w x y z
        A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
        vbar
        ebar
        edge_1_1 edge_2_1 edge_3_1 edge_4_1 edge_5_1 edge_6_1 edge_7_1 edge_8_1 edge_9_1 edge_10_1 edge_11_1 edge_12_1 edge_13_1 edge_14_1 edge_15_1 edge_16_1 edge_17_1 edge_18_1 edge_19_1 edge_20_1
        hedge0 hedge1 hedge2 hedge3 hedge4 hedge5 hedge6 hedge7 hedge8 hedge9 hedge10 hedge11 hedge12 hedge13 hedge14 hedge15 hedge16 hedge17 hedge18 hedge19 hedge20
        hedge_0 hedge_1 hedge_2 hedge_3 hedge_4 hedge_5 hedge_6 hedge_7 hedge_8 hedge_9 hedge_10 hedge_11 hedge_12 hedge_13 hedge_14 hedge_15 hedge_16 hedge_17 hedge_18 hedge_19 hedge_20 mul
    );

    pub fn test_initialize() {
        initialize();
        let _ = TS.p;
    }

    use core::panic;

    use crate::{
        network::{library::panicing::ErroringLibrary, parsing::NetworkParse},
        shadowing::symbolica_utils::TypstSettings,
        structure::{
            ToSymbolic, permuted,
            representation::{Euclidean, Lorentz, Minkowski, RepName, initialize},
            slot::IsAbstractSlot,
        },
        tensors::symbolic::SymbolicTensor,
    };

    use super::*;
    use insta::assert_snapshot;
    use library::DummyLibrary;
    use symbolica::{atom::AtomCore, parse, parse_lit, printer::PrintOptions, symbol};

    #[test]
    fn parse_scalar() {
        let expr = parse!("1");

        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    fn parse_pow() {
        let _sqrt = symbol!("sqrt_scalar", tag = SPENSO_TAG.tag);
        let expr = parse!("ee^3*sqrt_scalar((m+d(mink(4,1))^2)^(-5))");

        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        println!("{}", net.snapshot_dot());
        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    fn parse_ratio() {
        test_initialize();
        let expr = parse_lit!(
            (P(1, mink(4, 1)) * P(2, mink(4, 1))) / f(mul(P(3, mink(4, 1)) * P(4, mink(4, 1))))
        );
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_snapshot!(net.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "∏"];
          1	 [label = "S:(f(mul(P(3,mink(4,1))*P(4,mink(4,1)))))^(-1)"];
          2	 [label = "T:P(1,mink(4,1))"];
          3	 [label = "T:P(2,mink(4,1))"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
          2:6:s	-> 3:8:s	 [id=1 dir=none  color="red:blue;0.5" label="mink4|1"];
          2:5:s	-> 0:2:s	 [id=2  color="red:blue;0.5"];
          1:4:s	-> 0:3:s	 [id=3  color="red:blue;0.5"];
          3:7:s	-> 0:1:s	 [id=4  color="red:blue;0.5"];
        }
        "#);

        assert_eq!(net.simple_execute(), expr);
    }
    #[test]
    fn parse_div() {
        let expr = parse!(
            "c*a/bracket(d(mink(4,1))* b(mink(4,1)))(bracket(d(mink(4,1))* b(mink(4,1)))^-3)(bracket(d(mink(4,1))* b(mink(4,1)))^-2)"
        );
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_snapshot!(net.simple_execute().to_bare_ordered_string(), @"(b(mink(4,1)))^(-6)*(d(mink(4,1)))^(-6)*a*c");

        let expr = parse_lit!(st(Q(1, mink(4, 1)) * Q(2, mink(4, 1))) ^ -1);
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        println!("{}", net.snapshot_dot());
        assert_snapshot!(net.simple_execute().to_bare_ordered_string(), @"(st(Q(1,mink(4,1))*Q(2,mink(4,1))))^(-1)");
    }

    #[test]
    fn parse_scalar_tensor() {
        test_initialize();
        let expr = parse!("(
            (
                  -1*gammalooprs::{}::mUV^2+gammalooprs::{}::Q(6,spenso::mink(4,gammalooprs::{}::uv_mink_1337))
                    *gammalooprs::{}::Q(7,spenso::mink(4,gammalooprs::{}::uv_mink_1337))
                 )
            )*2");
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        println!("{}", net.dot_pretty());
        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    fn parse_val() {
        test_initialize();
        let expr = parse_lit!(
            ((Q(5, cind(0)))
                ^ 2 + (Q(5, cind(1)))
                ^ 2 * -1 + (Q(5, cind(2)))
                ^ 2 * -1 + (Q(5, cind(3)))
                ^ 2 * -1)
                ^ (-1)
                    * ((Q(6, cind(0)))
                        ^ 2 + (Q(6, cind(1)))
                        ^ 2 * -1 + (Q(6, cind(2)))
                        ^ 2 * -1 + (Q(6, cind(3)))
                        ^ 2 * -1)
                ^ (-1) * 1𝑖 / 3
                ^ 3 * Q(5, mink(4, edge_5_1))
                    * Q(6, mink(4, edge_6_1))
                    * u(1, bis(4, hedge_1))
                    * vbar(2, bis(4, hedge_2))
                    * ebar(0, mink(4, hedge_0))
                    * ebar(3, mink(4, hedge_3))
                    * ebar(4, mink(4, hedge_4))
                    * g(cof(3, hedge_1), dind(cof(3, hedge_2)))
                    * g(cof(3, hedge_1), dind(cof(3, hedge_2)))
                    * gamma(bis(4, hedge_2), bis(4, hedge_8), mink(4, hedge_0))
                    * gamma(bis(4, hedge_5), bis(4, hedge_1), mink(4, hedge_3))
                    * gamma(bis(4, hedge_6), bis(4, hedge_5), mink(4, edge_5_1))
                    * gamma(bis(4, hedge_7), bis(4, hedge_6), mink(4, hedge_4))
                    * gamma(bis(4, hedge_8), bis(4, hedge_7), mink(4, edge_6_1))
        );

        println!(
            "{}",
            expr.printer(PrintOptions {
                custom_print_mode: Some(("spenso", 0)),
                color_builtin_symbols: true,
                square_brackets_for_function: true,
                ..Default::default()
            })
        );
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_snapshot!(net.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "∏"];
          1	 [label = "S:((Q(5,cind(0)))^2+(Q(5,cind(1)))^2*-1+(Q(5,cind(2)))^2*-1+(Q(5,cind(3)))^2*-1)^(-1)*((Q(6,cind(0)))^2+(Q(6,cind(1)))^2*-1+(Q(6,cind(2)))^2*-1+(Q(6,cind(3)))^2*-1)^(-1)*(g(cof(3,hedge_1),dind(cof(3,hedge_2))))^2*1𝑖/27*u(1,bis(4,hedge_1))*vbar(2,bis(4,hedge_2))"];
          2	 [label = "T:Q(5,mink(4,edge_5_1))"];
          3	 [label = "T:Q(6,mink(4,edge_6_1))"];
          4	 [label = "T:ebar(0,mink(4,hedge_0))"];
          5	 [label = "T:ebar(3,mink(4,hedge_3))"];
          6	 [label = "T:ebar(4,mink(4,hedge_4))"];
          7	 [label = "T:gamma(bis(4,hedge_2),bis(4,hedge_8),mink(4,hedge_0))"];
          8	 [label = "T:gamma(bis(4,hedge_5),bis(4,hedge_1),mink(4,hedge_3))"];
          9	 [label = "T:gamma(bis(4,hedge_6),bis(4,hedge_5),mink(4,edge_5_1))"];
          10	 [label = "T:gamma(bis(4,hedge_7),bis(4,hedge_6),mink(4,hedge_4))"];
          11	 [label = "T:gamma(bis(4,hedge_8),bis(4,hedge_7),mink(4,edge_6_1))"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
          3:16:s	-> 11:32:s	 [id=1 dir=none  color="red:blue;0.5" label="mink4|edge_6_1"];
          6:22:s	-> 10:30:s	 [id=2 dir=none  color="red:blue;0.5" label="mink4|hedge_4"];
          2:14:s	-> 9:28:s	 [id=3 dir=none  color="red:blue;0.5" label="mink4|edge_5_1"];
          5:20:s	-> 8:26:s	 [id=4 dir=none  color="red:blue;0.5" label="mink4|hedge_3"];
          4:18:s	-> 7:24:s	 [id=5 dir=none  color="red:blue;0.5" label="mink4|hedge_0"];
          9:27:s	-> 0:3:s	 [id=6  color="red:blue;0.5"];
          7:23:s	-> 0:5:s	 [id=7  color="red:blue;0.5"];
          6:21:s	-> 0:6:s	 [id=8  color="red:blue;0.5"];
          10:29:s	-> 0:2:s	 [id=9  color="red:blue;0.5"];
          8:25:s	-> 0:4:s	 [id=10  color="red:blue;0.5"];
          1:12:s	-> 0:11:s	 [id=11  color="red:blue;0.5"];
          2:13:s	-> 0:10:s	 [id=12  color="red:blue;0.5"];
          3:15:s	-> 0:9:s	 [id=13  color="red:blue;0.5"];
          4:17:s	-> 0:8:s	 [id=14  color="red:blue;0.5"];
          5:19:s	-> 0:7:s	 [id=15  color="red:blue;0.5"];
          11:31:s	-> 0:1:s	 [id=16  color="red:blue;0.5"];
        }
        "#);
        assert_eq!(net.simple_execute(), expr);
        let mut out = String::new();
        expr.typst_fmt(&mut out, &TypstSettings::lowering())
            .unwrap();
        println!("{}", out)
    }

    #[test]
    fn parse_scalar_tensors_step_by() {
        let expr = parse!("c*a*b(mink(4,1))*d(mink(4,2))*d(mink(4,1))*d(mink(4,2))");

        let mut net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        let lib = DummyLibrary::<_>::new();
        let fnlib = ErroringLibrary::<Symbol>::new();

        let mut netc = net.clone();
        assert_snapshot!(
            net.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "∏"];
          1	 [label = "S:a*c"];
          2	 [label = "T:b(mink(4,1))"];
          3	 [label = "T:d(mink(4,1))"];
          4	 [label = "^( 2 )"];
          5	 [label = "T:d(mink(4,2))"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
          5:13:s	-> 4:11:s	 [id=1  color="red:blue;0.5"];
          2:7:s	-> 3:9:s	 [id=2 dir=none  color="red:blue;0.5" label="mink4|1"];
          2:6:s	-> 0:3:s	 [id=3  color="red:blue;0.5"];
          1:5:s	-> 0:4:s	 [id=4  color="red:blue;0.5"];
          3:8:s	-> 0:2:s	 [id=5  color="red:blue;0.5"];
          5:14:s	-> 4:12:s	 [id=6 dir=none  color="red:blue;0.5" label="mink4|2"];
          4:10:s	-> 0:1:s	 [id=7  color="red:blue;0.5"];
        }
        "#
        );
        net.execute::<Steps<1>, ContractScalars, _, _, _>(&lib, &fnlib)
            .unwrap();
        assert_snapshot!(
            net.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "∏"];
          1	 [label = "S:a*c"];
          2	 [label = "T:b(mink(4,1))"];
          3	 [label = "T:d(mink(4,1))"];
          4	 [label = "S:(d(mink(4,2)))^2"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
          3:8:s	-> 0:2:s	 [id=1  color="red:blue;0.5"];
          1:5:s	-> 0:4:s	 [id=2  color="red:blue;0.5"];
          2:6:s	-> 0:3:s	 [id=3  color="red:blue;0.5"];
          2:7:s	-> 3:9:s	 [id=4 dir=none  color="red:blue;0.5" label="mink4|1"];
          4:10:s	-> 0:1:s	 [id=5  color="red:blue;0.5"];
        }
        "#
        );
        net.execute::<Steps<1>, SingleSmallestDegree<false>, _, _, _>(&lib, &fnlib)
            .unwrap();
        assert_snapshot!(
            net.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          3	 [label = "∏"];
          0	 [label = "S:(d(mink(4,2)))^2"];
          1	 [label = "S:a*c"];
          2	 [label = "T:b(mink(4,1))*d(mink(4,1))"];
          ext0	 [style=invis];
          3:0:s	-> ext0	 [id=0 color="red"];
          1:5:s	-> 3:4:s	 [id=1  color="red:blue;0.5"];
          2:6:s	-> 3:3:s	 [id=2  color="red:blue;0.5"];
          0:2:s	-> 3:1:s	 [id=3  color="red:blue;0.5"];
        }
        "#
        );
        net.execute::<Steps<1>, SingleSmallestDegree<false>, _, _, _>(&lib, &fnlib)
            .unwrap();
        assert_snapshot!(
            net.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          3	 [label = "∏"];
          0	 [label = "S:(d(mink(4,2)))^2"];
          1	 [label = "S:a*c"];
          2	 [label = "T:b(mink(4,1))*d(mink(4,1))"];
          ext0	 [style=invis];
          3:0:s	-> ext0	 [id=0 color="red"];
          0:2:s	-> 3:1:s	 [id=1  color="red:blue;0.5"];
          1:5:s	-> 3:4:s	 [id=2  color="red:blue;0.5"];
          2:6:s	-> 3:3:s	 [id=3  color="red:blue;0.5"];
        }
        "#
        );
        net.execute::<Steps<1>, ContractScalars, _, _, _>(&lib, &fnlib)
            .unwrap();
        assert_snapshot!(
            net.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "T:(d(mink(4,2)))^2*a*b(mink(4,1))*c*d(mink(4,1))"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
        }
        "#
        );
        netc.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();
        assert_snapshot!(
            netc.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "T:(d(mink(4,2)))^2*a*b(mink(4,1))*c*d(mink(4,1))"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
        }
        "#
        );
        if let ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) =
            net.result().unwrap()
        {
            if let ExecutionResult::Val(TensorOrScalarOrKey::Tensor {
                tensor: tensor2, ..
            }) = netc.result().unwrap()
            {
                assert_eq!(tensor2.expression, tensor.expression);
            } else {
                panic!("Not scalar")
            }
        } else {
            panic!("Not scalar")
        }
    }

    #[test]
    fn parse_scalar_expr() {
        let expr = parse!("(y+x(mink(4,1))*y(mink(4,1))) *(1+1+2*x*(3*sin(r))/t)");
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    fn parse_tensor_expr() {
        let tensor1 = ShadowedStructure::<AbstractIndex>::from_iter(
            [
                Minkowski {}.new_slot(4, 1).to_lib(),
                Euclidean {}.new_slot(4, 2).to_lib(),
                Minkowski {}.new_slot(symbol!("d"), symbol!("mu")).to_lib(),
            ],
            symbol!("T"),
            None,
        )
        .structure
        .to_symbolic(None)
        .unwrap();

        let tensor2 = ShadowedStructure::<AbstractIndex>::from_iter(
            [
                Minkowski {}.new_slot(4, 1).to_lib(),
                Lorentz {}.new_slot(7, 1).to_lib(),
                Lorentz {}.new_slot(3, 2).to_lib(),
            ],
            symbol!("TT"),
            None,
        )
        .structure
        .to_symbolic(None)
        .unwrap();

        let tensor3 = ShadowedStructure::<AbstractIndex>::from_iter(
            [
                Lorentz {}.dual().new_slot(7, 1).to_lib(),
                Euclidean {}.new_slot(4, 2).to_lib(),
                Minkowski {}.new_slot(symbol!("d"), symbol!("mu")).to_lib(),
            ],
            symbol!("TTT"),
            None,
        )
        .structure
        .to_symbolic(None)
        .unwrap();

        let tensor4 = ShadowedStructure::<AbstractIndex>::from_iter(
            [Lorentz {}.new_slot(3, 2).to_lib()],
            symbol!("L"),
            None,
        )
        .structure
        .to_symbolic(None)
        .unwrap();

        let tensor5 = ShadowedStructure::<AbstractIndex>::from_iter(
            [Lorentz {}.dual().new_slot(3, 2).to_lib()],
            symbol!("P"),
            None,
        )
        .structure
        .to_symbolic(None)
        .unwrap();

        let expr = (parse!("a*sin(x/2)") * tensor1 * tensor2 * tensor3 + tensor4) * tensor5;

        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }

    // -G^2*(-g(mink(4,5),mink(4,6))*Q(2,mink(4,7))+g(mink(4,5),mink(4,6))*Q(3,mink(4,7))+g(mink(4,5),mink(4,7))*Q(2,mink(4,6))+g(mink(4,5),mink(4,7))*Q(4,mink(4,6))-g(mink(4,6),mink(4,7))*Q(3,mink(4,5))-g(mink(4,6),mink(4,7))*Q(4,mink(4,5)))*id(mink(4,2),mink(4,5))*id(mink(4,3),mink(4,6))*id(euc(4,0),euc(4,5))*id(euc(4,1),euc(4,4))*g(mink(4,4),mink(4,7))*vbar(1,euc(4,1))*u(0,euc(4,0))*ebar(2,mink(4,2))*ebar(3,mink(4,3))*gamma(mink(4,4),euc(4,5),euc(4,4))

    #[test]
    fn parse_big_tensors() {
        initialize();
        let expr = parse!(
            "-G^2*(-g(mink(4,5),mink(4,6))*Q(2,mink(4,7))+g(mink(4,5),mink(4,6))*Q(3,mink(4,7))+g(mink(4,5),mink(4,7))*Q(2,mink(4,6))+g(mink(4,5),mink(4,7))*Q(4,mink(4,6))-g(mink(4,6),mink(4,7))*Q(3,mink(4,5))-g(mink(4,6),mink(4,7))*Q(4,mink(4,5)))*id(mink(4,2),mink(4,5))*id(mink(4,3),mink(4,6))*id(euc(4,0),euc(4,5))*id(euc(4,1),euc(4,4))*g(mink(4,4),mink(4,7))*vbar(1,euc(4,1))*u(0,euc(4,0))*ebar(2,mink(4,2))*ebar(3,mink(4,3))*gamma(euc(4,5),euc(4,4),mink(4,4))"
        );
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    fn equal_duals() {
        initialize();
        let expr = parse_lit!(
            ((Q(5, cind(0)))
                ^ 2 + (Q(5, cind(1)))
                ^ 2 * -1 + (Q(5, cind(2)))
                ^ 2 * -1 + (Q(5, cind(3)))
                ^ 2 * -1)
                ^ (-1)
                    * ((Q(6, cind(0)))
                        ^ 2 + (Q(6, cind(1)))
                        ^ 2 * -1 + (Q(6, cind(2)))
                        ^ 2 * -1 + (Q(6, cind(3)))
                        ^ 2 * -1)
                ^ (-1)
                    * ((Q(7, cind(0)))
                        ^ 2 + (Q(7, cind(1)))
                        ^ 2 * -1 + (Q(7, cind(2)))
                        ^ 2 * -1 + (Q(7, cind(3)))
                        ^ 2 * -1)
                ^ (-1)
                    * ((Q(8, cind(0)))
                        ^ 2 + (Q(8, cind(1)))
                        ^ 2 * -1 + (Q(8, cind(2)))
                        ^ 2 * -1 + (Q(8, cind(3)))
                        ^ 2 * -1)
                ^ (-1)
                    * ((Q(9, cind(0)))
                        ^ 2 + (Q(9, cind(1)))
                        ^ 2 * -1 + (Q(9, cind(2)))
                        ^ 2 * -1 + (Q(9, cind(3)))
                        ^ 2 * -1)
                ^ (-1) * -1𝑖 / 3 * UFO::GC_11
                ^ 2 * UFO::GC_1
                ^ 3 * UFO::Gamma(mink(4, hedge_13), bis(4, hedge_11), 2)
                    * Q(5, mink(4, edge_5_1))
                    * Q(6, mink(4, edge_6_1))
                    * Q(7, mink(4, edge_7_1))
                    * Q(8, mink(4, edge_8_1))
                    * u(1)
                    * vbar(2, bis(4, hedge_2))
                    * ebar(0, mink(4, hedge_0))
                    * ebar(3, mink(4, hedge_3))
                    * ebar(4, mink(4, hedge_4))
                    * g(cof(3, hedge_1), dind(cof(3, hedge_2)))
                    * g(mink(4, hedge_13), mink(4, hedge_14))
                    * gamma(bis(4, hedge_10), bis(4, hedge_9), mink(4, edge_8_1))
                    * gamma(bis(4, hedge_12), bis(4, hedge_11), mink(4, edge_7_1))
                    * gamma(bis(4, hedge_2), bis(4, hedge_10), mink(4, hedge_14))
                    * gamma(bis(4, hedge_5), bis(4, hedge_12), mink(4, hedge_3))
                    * gamma(bis(4, hedge_6), bis(4, hedge_5), mink(4, edge_5_1))
                    * gamma(bis(4, hedge_7), bis(4, hedge_6), mink(4, hedge_4))
                    * gamma(bis(4, hedge_8), bis(4, hedge_7), mink(4, edge_6_1))
                    * gamma(bis(4, hedge_9), bis(4, hedge_8), mink(4, hedge_0))
                    * t(coad(8, hedge_13), cof(2), dind(cof(3, hedge_11)))
                    * t(coad(8, hedge_13), cof(3, hedge_11), dind(cof(3, hedge_2)))
        );

        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    fn gammaloop_six_photon() {
        initialize();
        let expr = parse_lit!(
            -64 / 243 * ee
                ^ 6 * (MT * g(euc(4, hedge3), euc(4, hedge4))
                    + gamma(euc(4, hedge3), euc(4, hedge4), mink(4, edge_3_1))
                        * Q(3, mink(4, edge_3_1)))
                    * (MT * g(euc(4, hedge6), euc(4, hedge7))
                        + gamma(euc(4, hedge6), euc(4, hedge7), mink(4, edge_5_1))
                            * Q(5, mink(4, edge_5_1)))
                    * (MT * g(euc(4, hedge9), euc(4, hedge10))
                        + gamma(euc(4, hedge9), euc(4, hedge10), mink(4, edge_7_1))
                            * Q(7, mink(4, edge_7_1)))
                    * (MT * g(euc(4, hedge11), euc(4, hedge12))
                        + gamma(euc(4, hedge11), euc(4, hedge12), mink(4, edge_8_1))
                            * Q(8, mink(4, edge_8_1)))
                    * (MT * g(euc(4, hedge13), euc(4, hedge14))
                        + gamma(euc(4, hedge13), euc(4, hedge14), mink(4, edge_9_1))
                            * Q(9, mink(4, edge_9_1)))
                    * (MT * g(euc(4, hedge16), euc(4, hedge17))
                        + gamma(euc(4, hedge16), euc(4, hedge17), mink(4, edge_11_1))
                            * Q(11, mink(4, edge_11_1)))
                    * gamma(euc(4, hedge4), euc(4, hedge6), mink(4, hedge5))
                    * gamma(euc(4, hedge7), euc(4, hedge9), mink(4, hedge8))
                    * gamma(euc(4, hedge10), euc(4, hedge11), mink(4, hedge0))
                    * gamma(euc(4, hedge12), euc(4, hedge13), mink(4, hedge1))
                    * gamma(euc(4, hedge14), euc(4, hedge16), mink(4, hedge15))
                    * gamma(euc(4, hedge17), euc(4, hedge3), mink(4, hedge2))
                    * eps(0, mink(4, hedge0))
                    * eps(1, mink(4, hedge1))
                    * epsbar(2, mink(4, hedge2))
                    * epsbar(4, mink(4, hedge5))
                    * epsbar(6, mink(4, hedge8))
                    * epsbar(10, mink(4, hedge15))
        );
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    fn parse_neg_tensors() {
        initialize();
        let expr =
            parse!("-d(mink(4,6),mink(4,5))*Q(2,mink(4,7))+d(mink(4,6),mink(4,5))*Q(3,mink(4,7))");
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    fn many_sums() {
        initialize();
        let expr = parse_lit!(
            (P(4, mink(4, r_2)) + N(4, mink(4, r_2)))
                * (P(5, mink(4, r_3)) + N(5, mink(4, r_3)))
                * (A(mink(4, r_2), mink(4, r_3)) + B(mink(4, r_3), mink(4, r_2)))
        );
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    fn contract_problem() {
        initialize();

        let expr = parse_lit!(
            (-1 * Q(EMRID(0, 4), mink(4, l_20))
                * gamma(euc(4, l_3), euc(4, l_6), mink(4, l_0))
                * gamma(euc(4, l_5), euc(4, l_2), mink(4, l_20))
                * gamma(euc(4, l_6), euc(4, l_5), mink(4, l_1))
                + 2 * Q(EMRID(0, 4), mink(4, l_1)) * gamma(euc(4, l_3), euc(4, l_2), mink(4, l_0)))
                * 1𝑖
                * G
                ^ 2
        );

        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }

    #[test]
    // #[should_panic]
    fn parse_problem() {
        initialize();
        let expr = parse_lit!(
            1 / 18 * ee
                ^ 2 * G
                ^ 4 * ((mUV
                    ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                        + OSE(3, Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                    ^ (1 / 2)
                        + (mUV
                            ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                                + OSE(4, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                    ^ (1 / 2))
                ^ -1 * ((mUV
                    ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                        + OSE(3, Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                    ^ (1 / 2)
                        + (mUV
                            ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                                + OSE(5, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                    ^ (1 / 2))
                ^ -1 * (-Q3(3, mink(4, edge_4_1))
                    + (mUV
                        ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                            + OSE(4, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                    ^ (1 / 2) * sigma(4) * delta(cind(0), mink(4, edge_4_1)))
                    * (-Q3(3, mink(4, edge_5_1))
                        + (mUV
                            ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                                + OSE(5, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                        ^ (1 / 2) * sigma(5) * delta(cind(0), mink(4, edge_5_1)))
                    * (mUV
                        ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                            + OSE(3, Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                ^ (-1 / 2)
                    * (mUV
                        ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                            + OSE(4, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                ^ (-1 / 2)
                    * (mUV
                        ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                            + OSE(5, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                ^ (-1 / 2)
                    * g(mink(4, hedge_3), mink(4, hedge_4))
                    * gamma(bis(4, hedge(1)), bis(4, hedge(8)), mink(4, hedge_4))
                ^ 2 * gamma(bis(4, hedge(5)), bis(4, hedge(0)), mink(4, hedge_3))
                ^ 2 * gamma(bis(4, hedge(6)), bis(4, hedge(5)), mink(4, edge_4_1))
                    * gamma(bis(4, hedge(7)), bis(4, hedge(6)), mink(4, hedge_2))
                ^ 2 * gamma(bis(4, hedge(8)), bis(4, hedge(7)), mink(4, edge_5_1)),
            default_namespace = "spenso"
        );

        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }
    #[test]
    // #[should_panic]
    fn infinite_execution() {
        initialize();
        let _s = symbol!("dim");

        let _ = parse_lit!(
            (N(0, mink(dim, l(0))) * P(0, mink(dim, r(0)))
                + N(0, mink(dim, r(0))) * P(0, mink(dim, l(0))))
                * N(0, mink(dim, dummy_ss(0, 1)))
                * P(0, mink(dim, dummy_ss(0, 1)))
                + -1 * N(0, mink(dim, dummy_ss(0, 3)))
                    * N(0, mink(dim, dummy_ss(0, 4)))
                    * P(0, mink(dim, dummy_ss(0, 3)))
                    * P(0, mink(dim, dummy_ss(0, 4)))
                    * g(mink(dim, l(0)), mink(dim, r(0)))
        )
        .replace(parse_lit!(dim))
        .with(Atom::num(4));

        let expr = parse_lit!(
            -1𝑖 * G
                ^ 3 * (g(mink(4, l_6), mink(4, l_8)) * g(mink(4, l_7), mink(4, l_9))
                    - g(mink(4, l_6), mink(4, l_9)) * g(mink(4, l_8), mink(4, l_7)))
                    * g(mink(4, l_0), mink(4, l_6))
                    * g(mink(4, l_1), mink(4, l_7))
                    * g(mink(4, l_4), mink(4, l_8))
                    * g(mink(4, l_5), mink(4, l_9))
                    * g(bis(4, l_2), bis(4, l_5))
                    * g(bis(4, l_3), bis(4, l_6))
                    * gamma(bis(4, l_6), bis(4, l_5), mink(4, l_5))
        );

        let lib = DummyLibrary::<_>::new();
        let fnlib = ErroringLibrary::<Symbol>::new();
        let mut net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        let mut net_iter = net.clone();

        loop {
            let old = net_iter.clone();
            net_iter
                .execute::<Steps<1>, SingleSmallestDegree<false>, _, _, _>(&lib, &fnlib)
                .unwrap();
            net_iter
                .execute::<Steps<1>, ContractScalars, _, _, _>(&lib, &fnlib)
                .unwrap();
            if net_iter == old {
                break;
            }
        }
        assert_snapshot!(
            net_iter.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "T:(-1*g(mink(4,l_6),mink(4,l_9))*g(mink(4,l_7),mink(4,l_8))+g(mink(4,l_6),mink(4,l_8))*g(mink(4,l_7),mink(4,l_9)))*-1𝑖*G^3*g(bis(4,l_2),bis(4,l_5))*g(bis(4,l_3),bis(4,l_6))*g(mink(4,l_0),mink(4,l_6))*g(mink(4,l_1),mink(4,l_7))*g(mink(4,l_4),mink(4,l_8))*g(mink(4,l_5),mink(4,l_9))*gamma(bis(4,l_6),bis(4,l_5),mink(4,l_5))"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
          ext1	 [style=invis];
          0:2:s	-> ext1	 [id=1 dir=none color="red" label="mink4|l_1"];
          ext2	 [style=invis];
          0:3:s	-> ext2	 [id=2 dir=none color="red" label="mink4|l_0"];
          ext3	 [style=invis];
          0:1:s	-> ext3	 [id=3 dir=none color="red" label="mink4|l_4"];
        }
        "#
        );
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();
    }

    #[test]
    fn gammaloop_input() {
        initialize();
        let expr = parse_lit!(
            16 / 81 * ee
                ^ 4 * (MT * g(euc(4, hedge_3), euc(4, hedge_4))
                    + gamma(euc(4, hedge_3), euc(4, hedge_4), mink(4, edge_3_1))
                        * Q(3, mink(4, edge_3_1)))
                    * (MT * g(euc(4, hedge_6), euc(4, hedge_7))
                    + gamma(euc(4, hedge_6), euc(4, hedge_7), mink(4, edge_5_1))
                            * Q(5, mink(4, edge_5_1)))
                    * (MT * g(euc(4, hedge_10), euc(4, hedge_11))
                    + gamma(euc(4, hedge_10), euc(4, hedge_11), mink(4, edge_8_1))
                            * Q(8, mink(4, edge_8_1)))
                    // * (-1 / 8 * (-Q(1, cind(0)) + Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (-Q(1, cind(0)) + Q(2, cind(0)) + Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (-Q(1, cind(0))
                    //         + Q(2, cind(0))
                    //         + Q(4, cind(0))
                    //         + Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, 1, 1, 1, 1, 1, 1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (-Q(1, cind(0)) + Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (-Q(1, cind(0)) + Q(2, cind(0)) + Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(4, cind(0)) - Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, 1, 1, -1, 1, 1, 1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (-Q(1, cind(0)) + Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (-Q(1, cind(0))
                    //         + Q(2, cind(0))
                    //         + Q(4, cind(0))
                    //         + Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, -1, 1, 1, 1, 1, 1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (-Q(1, cind(0)) + Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(4, cind(0)) - Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, -1, 1, -1, 1, 1, 1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (Q(1, cind(0)) - Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (-Q(1, cind(0)) + Q(2, cind(0)) + Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (-Q(1, cind(0))
                    //         + Q(2, cind(0))
                    //         + Q(4, cind(0))
                    //         + Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, 1, 1, 1, 1, 1, -1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (Q(1, cind(0)) - Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (-Q(1, cind(0)) + Q(2, cind(0)) + Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(4, cind(0)) - Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, 1, 1, -1, 1, 1, -1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (Q(1, cind(0)) - Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (-Q(1, cind(0))
                    //         + Q(2, cind(0))
                    //         + Q(4, cind(0))
                    //         + Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, -1, 1, 1, 1, 1, -1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (Q(1, cind(0)) - Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(4, cind(0)) - Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, -1, 1, -1, 1, 1, -1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1)
                    * g(dind(cof(3, hedge_8)), cof(3, hedge_1))
                    * gamma(euc(4, hedge_1), euc(4, hedge_10), mink(4, hedge_9))
                    * gamma(euc(4, hedge_4), euc(4, hedge_6), mink(4, hedge_5))
                    * gamma(euc(4, hedge_7), euc(4, hedge_8), mink(4, hedge_0))
                    * gamma(euc(4, hedge_11), euc(4, hedge_3), mink(4, hedge_2))
                    * ubar(6, euc(4, hedge_8))
                    * u(1, euc(4, hedge_1))
                    * eps(0, mink(4, hedge_0))
                    * eps(2, mink(4, hedge_2))
                    * eps(4, mink(4, hedge_5))
                    * eps(7, mink(4, hedge_9))
        );

        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(net.simple_execute(), expr);
    }
    #[test]
    fn wrapping() {
        initialize();
        let expr = parse_lit!(A * g(euc(4, hedge_3), euc(4, hedge_5)));

        let expr2 = parse_lit!(B * gg(euc(4, hedge_3), euc(4, hedge_4)));
        let expr3 = parse_lit!(C * ggg(euc(4, hedge_5), euc(4, hedge_4)));

        let lib = DummyLibrary::<SymbolicTensor>::new();
        let fnlib = ErroringLibrary::<Symbol>::new();
        let settings = &ParseSettings::default();
        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(settings)
            .unwrap();
        let net2 = expr2
            .parse_to_symbolic_net::<AbstractIndex>(settings)
            .unwrap();
        let net3 = expr3
            .parse_to_symbolic_net::<AbstractIndex>(settings)
            .unwrap();

        let mut acc = Network::one();
        println!("{}", expr);

        acc *= net;
        acc *= net2;
        acc *= net3;
        println!("{}", acc.snapshot_dot());

        acc.merge_ops();

        println!("{}", acc.snapshot_dot());

        // acc.execute::<Sequential, SmallestDegree, _, _,_>(&lib,&fnlib)
        //     .unwrap();

        acc.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();
        println!("{}", acc.snapshot_dot());
        let obt: Atom = acc.result_scalar().unwrap().into();
        assert_eq!(obt, expr * expr2 * expr3)
    }
    #[test]
    fn scalar_mult() {
        initialize();
        let expr = parse_lit!(A * B * C * g(euc(4, hedge_3)));

        let expr2 = parse_lit!(B * gg(euc(4, hedge_3)));
        let expr3 = parse_lit!(C * ggg(euc(4, hedge_5)) * g(euc(4, hedge_4)));
        let expr4 = parse_lit!(A * B * ggg(euc(4, hedge_5)) * g(euc(4, hedge_4)));
        let fnlib = ErroringLibrary::<Symbol>::new();
        let lib = DummyLibrary::<SymbolicTensor>::new();

        for ex in [expr, expr2, expr3, expr4] {
            println!("Expr:{ex}");
            let mut acc = ex
                .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap();
            // acc *= net3;
            println!("{}", acc.snapshot_dot());

            acc.execute::<Steps<1>, ContractScalars, _, _, _>(&lib, &fnlib)
                .unwrap();

            println!("{}", acc.snapshot_dot());
            acc.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
                .unwrap();

            println!("{}", acc.snapshot_dot());

            if let ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) =
                acc.result().unwrap()
            {
                // println!("YaY:{}", (&expr - &tensor.expression).expand());
                assert_eq!(ex, tensor.expression);
            } else {
                panic!("Not tensor")
            }
        }
    }

    #[test]
    fn dot_derivative() {
        initialize();

        let l = symbol!("lambda";Scalar);
        let expr = parse_lit!(
            // (2 * dot(Q3(0), Q3(2), mink(4)) + dot(Q3(0), Q3(0), mink(4)))
            dot(mink(4), Q(1), Q(3)) * dot(mink(4), Q(3), Q(1) + Q(2))
        );

        // 2*dot(Q3(0),Q3(2),mink(4)) is Scalar vs dot(Q3(0),Q3(0),mink(4)

        dbg!(<permuted::PermutedStructure<
            NamedStructure<symbolica::atom::Symbol, Vec<symbolica::atom::Atom>>,
        > as TryFrom<AtomView>>::Error::ParsingError(
            "a".into()
        ));

        let e = expr
            .replace(parse_lit!(Q(1)))
            .with(parse_lit!(Q(1)) * l)
            .expand();
        println!("{}", e);

        let net = e
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_snapshot!(net.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "∑"];
          1	 [label = "∏"];
          2	 [label = "S:lambda^2"];
          3	 [label = "^( 2 )"];
          4	 [label = "∏"];
          5	 [label = "T:Q(1,mink(4,1))"];
          6	 [label = "T:Q(3,mink(4,1))"];
          7	 [label = "∏"];
          8	 [label = "S:lambda"];
          9	 [label = "∏"];
          10	 [label = "T:Q(1,mink(4,1))"];
          11	 [label = "T:Q(3,mink(4,1))"];
          12	 [label = "∏"];
          13	 [label = "T:Q(2,mink(4,1))"];
          14	 [label = "T:Q(3,mink(4,1))"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
          7:16:s	-> 0:1:s	 [id=1  color="red:blue;0.5"];
          1:3:s	-> 0:2:s	 [id=2  color="red:blue;0.5"];
          6:14:s	-> 4:10:s	 [id=3  color="red:blue;0.5"];
          3:7:s	-> 1:4:s	 [id=4  color="red:blue;0.5"];
          2:6:s	-> 1:5:s	 [id=5  color="red:blue;0.5"];
          5:12:s	-> 4:11:s	 [id=6  color="red:blue;0.5"];
          4:9:s	-> 3:8:s	 [id=7  color="red:blue;0.5"];
          5:13:s	-> 6:15:s	 [id=8 dir=none  color="red:blue;0.5" label="mink4|1"];
          13:31:s	-> 12:30:s	 [id=9  color="red:blue;0.5"];
          12:28:s	-> 7:17:s	 [id=10  color="red:blue;0.5"];
          9:21:s	-> 7:18:s	 [id=11  color="red:blue;0.5"];
          8:20:s	-> 7:19:s	 [id=12  color="red:blue;0.5"];
          10:25:s	-> 11:27:s	 [id=13 dir=none  color="red:blue;0.5" label="mink4|1"];
          11:26:s	-> 9:22:s	 [id=14  color="red:blue;0.5"];
          10:24:s	-> 9:23:s	 [id=15  color="red:blue;0.5"];
          13:32:s	-> 14:34:s	 [id=16 dir=none  color="red:blue;0.5" label="mink4|1"];
          14:33:s	-> 12:29:s	 [id=17  color="red:blue;0.5"];
        }
        "#);

        assert_snapshot!(net.simple_execute().to_bare_ordered_string(),@"(Q(1,mink(4,1)))^2*(Q(3,mink(4,1)))^2*lambda^2+(Q(3,mink(4,1)))^2*Q(1,mink(4,1))*Q(2,mink(4,1))*lambda");
    }
}
