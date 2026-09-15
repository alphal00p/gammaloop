use eyre::eyre;
use std::collections::HashMap;
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol},
    function,
    printer::PrintState,
    symbol,
};

use colored::Colorize;

use crate::{
    network::{
        library::{FunctionLibrary, FunctionLibraryError},
        tags::SPENSO_TAG,
    },
    structure::{HasStructure, TensorStructure},
    tensors::{
        data::StorageTensor,
        parametric::{ParamOrConcrete, ParamTensor, to_param::ToParam},
    },
};

pub struct Inbuilts {
    pub conj: Symbol,
}

impl Inbuilts {
    pub fn conj<'a, A: Into<AtomOrView<'a>>>(&self, atom: A) -> Atom {
        let a = atom.into();
        function!(self.conj, a.as_view())
    }
}
crate::symbolica_init_lazy_static! {
    pub static INBUILTS, INBUILTS_INNER: Inbuilts = || Inbuilts {
        conj: symbol!(
            "spenso::conj",
            tag = SPENSO_TAG.broadcast,
            norm = |view, out| {
                if let AtomView::Fun(dind1) = view
                    && dind1.get_nargs() == 1
                {
                    let arg = dind1.iter().next().unwrap();
                    if let AtomView::Fun(arg) = arg
                        && arg.get_nargs() == 1
                        && arg.get_symbol() == symbol!("spenso::conj")
                    {
                        **out = arg.iter().next().unwrap().to_owned();
                    }
                }
            },
            print = |a, opt, _state| {
                if opt.color_builtin_symbols {
                    let mut fmt = "conj".blue().to_string();
                    if let AtomView::Fun(f) = a {
                        fmt.push('(');
                        let n_args = f.get_nargs();
                        for (i, a) in f.iter().enumerate() {
                            a.format(&mut fmt, opt, PrintState::new()).unwrap();
                            if i < n_args - 1 {
                                fmt.push(',');
                            }
                        }
                        fmt.push(')');
                    }

                    Some(fmt)
                } else {
                    None
                }
            }
        ),
    };
}

pub struct SymbolLib<T, Missing> {
    #[allow(clippy::type_complexity)]
    pub functions:
        HashMap<Symbol, Box<dyn Fn(T) -> Result<T, FunctionLibraryError<Symbol>> + Send + Sync>>,
    #[allow(clippy::type_complexity)]
    pub scalar_functions: HashMap<
        Symbol,
        Box<dyn Fn(Atom) -> Result<Atom, FunctionLibraryError<Symbol>> + Send + Sync>,
    >,
    pub _missing: Missing,
}

impl<T, Missing> SymbolLib<T, Missing> {
    pub fn insert<F>(&mut self, key: Symbol, func: F)
    where
        F: Fn(T) -> T + Send + Sync + 'static,
    {
        self.insert_fallible(key, move |tensor| Ok(func(tensor)));
    }

    pub fn insert_fallible<F>(&mut self, key: Symbol, func: F)
    where
        F: Fn(T) -> Result<T, FunctionLibraryError<Symbol>> + Send + Sync + 'static,
    {
        self.functions.insert(key, Box::new(func));
    }

    pub fn insert_scalar_fallible<F>(&mut self, key: Symbol, func: F)
    where
        F: Fn(Atom) -> Result<Atom, FunctionLibraryError<Symbol>> + Send + Sync + 'static,
    {
        self.scalar_functions.insert(key, Box::new(func));
    }
}

pub struct Panic;
impl Panic {
    pub fn new_lib<T>() -> SymbolLib<T, Self> {
        SymbolLib {
            functions: HashMap::new(),
            scalar_functions: HashMap::new(),
            _missing: Self,
        }
    }
}

impl<S: TensorStructure> FunctionLibrary<ParamTensor<S>, Atom>
    for SymbolLib<ParamTensor<S>, Panic>
{
    type Key = Symbol;

    fn apply(
        &self,
        key: &Self::Key,
        tensor: ParamTensor<S>,
    ) -> Result<ParamTensor<S>, FunctionLibraryError<Symbol>> {
        if let Some(func) = self.functions.get(key) {
            func(tensor)
        } else {
            Err(FunctionLibraryError::NotFound(*key))
        }
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: Atom,
    ) -> eyre::Result<Atom, FunctionLibraryError<Self::Key>> {
        if let Some(function) = self.scalar_functions.get(key) {
            function(scalar)
        } else {
            Ok(function!(*key, scalar))
        }
    }
}

pub struct Wrap;
impl Wrap {
    pub fn new_lib<T>() -> SymbolLib<T, Self> {
        SymbolLib {
            functions: HashMap::new(),
            scalar_functions: HashMap::new(),
            _missing: Self,
        }
    }
}
impl<S: TensorStructure + Clone> FunctionLibrary<ParamTensor<S>, Atom>
    for SymbolLib<ParamTensor<S>, Wrap>
{
    type Key = Symbol;

    fn apply(
        &self,
        key: &Self::Key,
        tensor: ParamTensor<S>,
    ) -> Result<ParamTensor<S>, FunctionLibraryError<Symbol>> {
        if let Some(func) = self.functions.get(key) {
            func(tensor)
        } else {
            Ok(tensor.map_data_self(|a| function!(*key, a)))
        }
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: Atom,
    ) -> eyre::Result<Atom, FunctionLibraryError<Self::Key>> {
        if let Some(function) = self.scalar_functions.get(key) {
            function(scalar)
        } else {
            Ok(function!(*key, scalar))
        }
    }
}

impl<S: TensorStructure + Clone, C: ToParam + HasStructure<Structure = S>>
    FunctionLibrary<ParamOrConcrete<C, S>, Atom> for SymbolLib<C, Wrap>
{
    type Key = Symbol;

    fn apply(
        &self,
        key: &Self::Key,
        tensor: ParamOrConcrete<C, S>,
    ) -> Result<ParamOrConcrete<C, S>, FunctionLibraryError<Symbol>> {
        Ok(match tensor {
            ParamOrConcrete::Concrete(c) => {
                if let Some(func) = self.functions.get(key) {
                    ParamOrConcrete::Concrete(func(c)?)
                } else {
                    ParamOrConcrete::Param(c.to_param().map_data_self(|a| function!(*key, a)))
                }
            }
            ParamOrConcrete::Param(p) => {
                ParamOrConcrete::Param(p.map_data_self(|a| function!(*key, a)))
            }
        })
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: Atom,
    ) -> eyre::Result<Atom, FunctionLibraryError<Self::Key>> {
        if let Some(function) = self.scalar_functions.get(key) {
            function(scalar)
        } else {
            Ok(function!(*key, scalar))
        }
    }
}

pub struct PanicMissingConcrete;
impl PanicMissingConcrete {
    pub fn new_lib<T>() -> SymbolLib<T, Self> {
        SymbolLib {
            functions: HashMap::new(),
            scalar_functions: HashMap::new(),
            _missing: Self,
        }
    }
}

impl<S: TensorStructure + Clone, C: ToParam + HasStructure<Structure = S>>
    FunctionLibrary<ParamOrConcrete<C, S>, Atom> for SymbolLib<C, PanicMissingConcrete>
{
    type Key = Symbol;

    fn apply(
        &self,
        key: &Self::Key,
        tensor: ParamOrConcrete<C, S>,
    ) -> Result<ParamOrConcrete<C, S>, FunctionLibraryError<Symbol>> {
        match tensor {
            ParamOrConcrete::Concrete(c) => {
                if let Some(func) = self.functions.get(key) {
                    Ok(ParamOrConcrete::Concrete(func(c)?))
                } else {
                    Err(FunctionLibraryError::Other(eyre!(
                        "no concrete callback registered for broadcast function `{key}`"
                    )))
                }
            }
            ParamOrConcrete::Param(p) => Ok(ParamOrConcrete::Param(
                p.map_data_self(|a| function!(*key, a)),
            )),
        }
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: Atom,
    ) -> eyre::Result<Atom, FunctionLibraryError<Self::Key>> {
        if let Some(function) = self.scalar_functions.get(key) {
            function(scalar)
        } else if symbolica::domains::float::Complex::<f64>::try_from(scalar.as_view()).is_ok() {
            Err(FunctionLibraryError::Other(eyre!(
                "no concrete callback registered for broadcast function `{key}`"
            )))
        } else {
            Ok(function!(*key, scalar))
        }
    }
}

impl<S: TensorStructure + Clone, C: ToParam + HasStructure<Structure = S>>
    FunctionLibrary<ParamOrConcrete<C, S>, Atom> for SymbolLib<C, Panic>
{
    type Key = Symbol;

    fn apply(
        &self,
        key: &Self::Key,
        tensor: ParamOrConcrete<C, S>,
    ) -> Result<ParamOrConcrete<C, S>, FunctionLibraryError<Symbol>> {
        if let Some(func) = self.functions.get(key) {
            if let ParamOrConcrete::Concrete(c) = tensor {
                Ok(ParamOrConcrete::Concrete(func(c)?))
            } else {
                Err(FunctionLibraryError::Other(eyre!(
                    "Cannot map parametric tensor"
                )))
            }
        } else {
            Err(FunctionLibraryError::Other(eyre!(
                "no concrete callback registered for broadcast function `{key}`"
            )))
        }
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: Atom,
    ) -> eyre::Result<Atom, FunctionLibraryError<Self::Key>> {
        if let Some(function) = self.scalar_functions.get(key) {
            function(scalar)
        } else {
            Ok(function!(*key, scalar))
        }
    }
}

impl<S: TensorStructure + Clone, C: ToParam + HasStructure<Structure = S>>
    FunctionLibrary<ParamOrConcrete<C, S>, Atom> for SymbolLib<ParamOrConcrete<C, S>, Wrap>
{
    type Key = Symbol;

    fn apply(
        &self,
        key: &Self::Key,
        tensor: ParamOrConcrete<C, S>,
    ) -> Result<ParamOrConcrete<C, S>, FunctionLibraryError<Symbol>> {
        if let Some(func) = self.functions.get(key) {
            func(tensor)
        } else {
            Ok(ParamOrConcrete::Param(
                match tensor {
                    ParamOrConcrete::Concrete(c) => c.to_param(),
                    ParamOrConcrete::Param(p) => p,
                }
                .map_data_self(|a| function!(*key, a)),
            ))
        }
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: Atom,
    ) -> eyre::Result<Atom, FunctionLibraryError<Self::Key>> {
        if let Some(function) = self.scalar_functions.get(key) {
            function(scalar)
        } else {
            Ok(function!(*key, scalar))
        }
    }
}

impl<S: TensorStructure + Clone, C> FunctionLibrary<ParamOrConcrete<C, S>, Atom>
    for SymbolLib<ParamOrConcrete<C, S>, PanicMissingConcrete>
{
    type Key = Symbol;

    fn apply(
        &self,
        key: &Self::Key,
        tensor: ParamOrConcrete<C, S>,
    ) -> Result<ParamOrConcrete<C, S>, FunctionLibraryError<Symbol>> {
        if let Some(func) = self.functions.get(key) {
            func(tensor)
        } else if let ParamOrConcrete::Param(p) = tensor {
            Ok(ParamOrConcrete::Param(
                p.map_data_self(|a| function!(*key, a)),
            ))
        } else {
            Err(FunctionLibraryError::NotFound(*key))
        }
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: Atom,
    ) -> eyre::Result<Atom, FunctionLibraryError<Self::Key>> {
        if let Some(function) = self.scalar_functions.get(key) {
            function(scalar)
        } else if symbolica::domains::float::Complex::<f64>::try_from(scalar.as_view()).is_ok() {
            Err(FunctionLibraryError::Other(eyre!(
                "no concrete callback registered for broadcast function `{key}`"
            )))
        } else {
            Ok(function!(*key, scalar))
        }
    }
}

impl<S: TensorStructure + Clone, C: ToParam + HasStructure<Structure = S>>
    FunctionLibrary<ParamOrConcrete<C, S>, Atom> for SymbolLib<ParamOrConcrete<C, S>, Panic>
{
    type Key = Symbol;
    fn apply(
        &self,
        key: &Self::Key,
        tensor: ParamOrConcrete<C, S>,
    ) -> Result<ParamOrConcrete<C, S>, FunctionLibraryError<Symbol>> {
        if let Some(func) = self.functions.get(key) {
            func(tensor)
        } else {
            Err(FunctionLibraryError::NotFound(*key))
        }
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: Atom,
    ) -> eyre::Result<Atom, FunctionLibraryError<Self::Key>> {
        if let Some(function) = self.scalar_functions.get(key) {
            function(scalar)
        } else {
            Ok(function!(*key, scalar))
        }
    }
}

#[test]
fn conj_construction() {
    // let a=  symbol!("spenso::conj", tag = SPENSO_TAG.tag);
}

#[test]
fn infallible_function_adapter_preserves_results() {
    let key = symbol!("spenso::test_infallible_function_adapter");
    let mut library = Panic::new_lib();
    library.insert(key, |value: i32| value + 1);

    assert_eq!(library.functions[&key](1).unwrap(), 2);
}

#[test]
fn fallible_function_adapter_preserves_errors() {
    let key = symbol!("spenso::test_fallible_function_adapter");
    let mut library = Panic::new_lib();
    library.insert_fallible(key, |_value: i32| {
        Err(FunctionLibraryError::Other(eyre!("callback failed")))
    });

    assert_eq!(
        library.functions[&key](1).unwrap_err().to_string(),
        "callback failed"
    );
}
