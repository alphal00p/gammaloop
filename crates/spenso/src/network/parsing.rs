use symbolica::atom::{AtomCore, AtomOrView, FunctionBuilder, Symbol};
use symbolica::domains::rational::Rational;
use symbolica::printer::PrintState;
use symbolica::{symbol, tag};

use super::*;

use library::Library;

use crate::network::library::DummyLibrary;
use crate::network::library::panicing::ErroringLibrary;
use crate::shadowing::symbolica_utils::SpensoPrintSettings;

use crate::network::library::symbolic::ETS;
use crate::structure::abstract_index::AIND_SYMBOLS;
use crate::structure::representation::Representation;
// use crate::shadowing::Concretize;
use crate::structure::slot::{DummyAind, ParseableAind, Slot, SlotError};
use crate::structure::{
    NamedStructure, OrderedStructure, PermutedStructure, StructureContract, StructureError,
    TensorShell,
};
use crate::tensors::parametric::ParamTensor;

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
    pub tensor: String,
    pub index: String,
    pub representation: String,
    pub i_: Symbol,
    pub dot: Symbol,
    pub rep_: Symbol,
    pub self_dual: String,
    pub self_dual_: Symbol,
    pub dualizable: String,
    pub dualizable_: Symbol,
}

impl SpensoTags {
    pub fn self_dual_<'a, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        let mut f = FunctionBuilder::new(self.self_dual_);
        for a in args.into_iter() {
            f = f.add_arg(a);
        }
        f.finish()
    }

    pub fn dualizable_<'a, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        let mut f = FunctionBuilder::new(self.dualizable_);
        for a in args.into_iter() {
            f = f.add_arg(a);
        }
        f.finish()
    }

    pub fn dualizable_dual_<'a, A: Into<AtomOrView<'a>>>(
        &self,
        args: impl IntoIterator<Item = A>,
    ) -> Atom {
        AIND_SYMBOLS.dual(self.dualizable_(args))
    }
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
                with_dim,..
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
            if with_dim {a.format(&mut s, opt, PrintState::new()).unwrap();
                s.push('.');
            }
            c.format(&mut s, opt,PrintState::new()).unwrap();
            if parens {
                s.push(')');
            }
            Some(s)

        },
        _=>None
    }



    }),

    tensor: tag!("tensor"),
    index: tag!("index"),
    self_dual: tag!("self_dual"),
    dualizable: tag!("dualizable"),
    representation: tag!("representation"),
    i_: symbol!("i_", tag = &tag!("index")),
    rep_: symbol!("rep_", tag = &tag!("representation")),

    self_dual_: symbol!(
        "self_dual_",
        tags = [&tag!("self_dual"), &tag!("representation")]
    ),
    dualizable_: symbol!(
        "dualizable_",
        tags = [&tag!("dualizable"), &tag!("representation")]
    ),
});

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

pub trait Parse: Sized {
    fn parse(view: AtomView) -> Result<PermutedStructure<Self>, StructureError>;
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
    fn as_leaf<S>(value: AtomView<'a>) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let s: Result<PermutedStructure<S>, _> = S::parse(value);

        println!("Looking at :{}", value);
        return if let Ok(s) = s {
            println!("Tensor");
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
        println!("Mul");
        if let Some(a) = settings.depth_limit
            && a <= state.depth
        {
            // println!("Mul leaf");
            return Self::as_leaf::<S>(value.as_view());
        }

        state.depth += 1;
        println!("{} for mul {}", state.depth, value.as_view());
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
        } else if symbol == SPENSO_TAG.pure_scalar {
            if value.get_nargs() != 1 {
                return Err(TensorNetworkError::TooManyArgsFunction(
                    value.as_view().to_plain_string(),
                ));
            }

            Ok(Self::from_scalar(value.iter().next().unwrap().try_into()?))
        } else if symbol == ETS.metric && settings.parse_inner_products && value.get_nargs() == 3 {
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
            let s: Result<PermutedStructure<S>, _> = S::parse(value.as_view());

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
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        if let Some(a) = settings.depth_limit
            && a <= state.depth
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
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        if let Some(a) = settings.depth_limit
            && a <= state.depth
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
