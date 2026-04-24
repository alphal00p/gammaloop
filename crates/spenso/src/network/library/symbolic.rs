use super::*;
use ahash::AHashMap;
use eyre::eyre;
use linnet::permutation::Permutation;
use std::sync::LazyLock;

use symbolica::printer::PrintState;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    symbol,
};

use crate::shadowing::symbolica_utils::SpensoPrintSettings;
use crate::{
    shadowing::symbolica_utils::{IntoArgs, IntoSymbol},
    structure::{
        HasName, IndexlessNamedStructure,
        abstract_index::AIND_SYMBOLS,
        named::{IdentityName, METRIC_NAME},
        permuted::{Perm, PermuteTensor},
        representation::{LibraryRep, RepName, initialize},
        slot::AbsInd,
    },
    tensors::parametric::{ConcreteOrParam, MixedTensor, ParamOrConcrete, ParamTensor},
};

pub type ExplicitKey<Aind> = IndexlessNamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind>;
pub type LibraryKey<Aind> = PermutedStructure<ExplicitKey<Aind>>;
impl<Aind> ExplicitKey<Aind> {
    pub fn from_structure<S: TensorStructure + HasName<Name: IntoSymbol, Args: IntoArgs>>(
        structure: &PermutedStructure<S>,
    ) -> Option<Self> {
        let rep_structure: Vec<_> = structure
            .structure
            .reps()
            .into_iter()
            .map(|r| r.to_lib())
            .collect();

        Some(
            IndexlessNamedStructure::from_iter(
                rep_structure,
                structure.structure.name()?.ref_into_symbol(),
                structure.structure.args().map(|a| a.args()),
            )
            .structure,
        )
    }
}

// pub struct DataStoreRefTensor<'a, T, S> {
//     data: &'a T,
//     structure: S,
// }

impl<S: TensorStructure + Clone> LibraryTensor for ParamTensor<S> {
    type WithIndices = ParamTensor<S::Indexed>;
    type Data = Atom;
    fn empty(key: S, zero: Atom) -> Self {
        ParamTensor::composite(<DataTensor<Atom, S> as LibraryTensor>::empty(key, zero))
    }

    fn from_dense(key: S, data: Vec<Self::Data>) -> Result<Self> {
        Ok(ParamTensor::composite(
            <DataTensor<Atom, S> as LibraryTensor>::from_dense(key, data)?,
        ))
    }

    fn from_sparse(
        key: S,
        data: impl IntoIterator<Item = (Vec<ConcreteIndex>, Self::Data)>,
        zero: Atom,
    ) -> Result<Self> {
        Ok(ParamTensor::composite(
            <DataTensor<Atom, S> as LibraryTensor>::from_sparse(key, data, zero)?,
        ))
    }

    fn with_indices(
        &self,
        indices: &[<<<Self::WithIndices as HasStructure>::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::WithIndices>, StructureError> {
        let new_tensor =
            <DataTensor<Atom, S> as LibraryTensor>::with_indices(&self.tensor, indices)?;
        Ok(PermutedStructure {
            structure: ParamTensor {
                tensor: new_tensor.structure,
                param_type: self.param_type,
            },
            rep_permutation: new_tensor.rep_permutation,
            index_permutation: new_tensor.index_permutation,
        })
    }
}

impl<D: Default + Clone, S: TensorStructure + Clone> LibraryTensor for MixedTensor<D, S> {
    type WithIndices = MixedTensor<D, S::Indexed>;
    type Data = ConcreteOrParam<RealOrComplex<D>>;

    fn empty(key: S, zero: Self::Data) -> Self {
        MixedTensor::Concrete(<RealOrComplexTensor<D, S> as LibraryTensor>::empty(
            key,
            match zero {
                ConcreteOrParam::Concrete(c) => c,
                ConcreteOrParam::Param(_) => RealOrComplex::Real(D::default()),
            },
        ))
    }

    fn from_dense(key: S, data: Vec<Self::Data>) -> Result<Self> {
        let data: Result<Vec<_>> = data
            .into_iter()
            .map(|v| match v {
                ConcreteOrParam::Concrete(c) => Ok(c),
                ConcreteOrParam::Param(p) => Err(eyre!("Only concrete data allowed, not {p}")),
            })
            .collect();

        Ok(MixedTensor::Concrete(
            <RealOrComplexTensor<D, S> as LibraryTensor>::from_dense(key, data?)?,
        ))
    }

    fn from_sparse(
        key: S,
        data: impl IntoIterator<Item = (Vec<ConcreteIndex>, Self::Data)>,
        zero: Self::Data,
    ) -> Result<Self> {
        let data: Result<Vec<_>> = data
            .into_iter()
            .map(|(i, v)| match v {
                ConcreteOrParam::Concrete(c) => Ok((i, c)),
                ConcreteOrParam::Param(p) => Err(eyre!("Only concrete data allowed, not {p}")),
            })
            .collect();

        Ok(MixedTensor::Concrete(
            <RealOrComplexTensor<D, S> as LibraryTensor>::from_sparse(
                key,
                data?,
                match zero {
                    ConcreteOrParam::Concrete(c) => c,
                    ConcreteOrParam::Param(_) => RealOrComplex::Real(D::default()),
                },
            )?,
        ))
    }

    fn with_indices(
        &self,
        indices: &[<<<Self::WithIndices as HasStructure>::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::WithIndices>, StructureError> {
        Ok(match self {
            ParamOrConcrete::Concrete(c) => {
                let strct = <RealOrComplexTensor<D, S> as LibraryTensor>::with_indices(c, indices)?;
                PermutedStructure {
                    structure: ParamOrConcrete::Concrete(strct.structure),
                    rep_permutation: strct.rep_permutation,
                    index_permutation: strct.index_permutation,
                }
            }
            ParamOrConcrete::Param(p) => {
                let strct = <ParamTensor<S> as LibraryTensor>::with_indices(p, indices)?;
                PermutedStructure {
                    structure: ParamOrConcrete::Param(strct.structure),
                    rep_permutation: strct.rep_permutation,
                    index_permutation: strct.index_permutation,
                }
            }
        })
    }
}

#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub struct GenericKey {
    global_name: Symbol,
    args: Option<Vec<Atom>>,
    reps: Vec<LibraryRep>,
}

impl GenericKey {
    pub fn new(global_name: Symbol, args: Option<Vec<Atom>>, reps: Vec<LibraryRep>) -> Self {
        Self {
            global_name,
            args,
            reps,
        }
    }
}

impl<Aind: AbsInd> From<ExplicitKey<Aind>> for GenericKey {
    fn from(key: ExplicitKey<Aind>) -> Self {
        Self {
            global_name: key.global_name.unwrap(),
            reps: key.reps().into_iter().map(|r| r.rep).collect(),
            args: key.additional_args,
        }
    }
}

#[allow(clippy::type_complexity)]
pub struct TensorLibrary<T: HasStructure<Structure = ExplicitKey<Aind>>, Aind> {
    explicit_dimension: AHashMap<ExplicitKey<Aind>, PermutedStructure<T>>,
    generic_dimension: AHashMap<GenericKey, fn(ExplicitKey<Aind>) -> T>,
}

pub struct ExplicitTensorSymbols {
    pub flat: Symbol,
    /// Print modes:
    /// - spenso: coloring builtin symbols colors the representation name, and tensor name
    ///   - 0: no dim,
    ///   - 1: no dim with parenthesis,
    ///   - 2: no dim no commas
    ///   - 3: no dim with parenthesis no commas
    ///   - 4: with dim,
    ///   - 5: with dim and parenthesis
    ///   - 6: with dim no commas
    ///   - 7: with dim and parenthesis no commas
    pub metric: Symbol,
}

pub static ETS: LazyLock<ExplicitTensorSymbols> = LazyLock::new(|| ExplicitTensorSymbols {
    flat: symbol!("♭";Symmetric;print = |a, opt| {

        match opt.custom_print_mode {
            Some(("spenso",i))=>{
                let SpensoPrintSettings{
                    parens,
                    commas,..
                } = SpensoPrintSettings::from(i);


                let AtomView::Fun(f)=a else {
                    return None;
                };
                if f.get_nargs()!=2 {
                    return None;
                }
                let mut argitem = f.iter();
                let a = argitem.next().unwrap();
                let b = argitem.next().unwrap();

                let mut out = if opt.color_builtin_symbols {
                    nu_ansi_term::Color::Magenta.paint("♭_").to_string()
                } else {
                    "♭_".to_string()
                };

                if parens{
                    out.push('(');
                }
                a.format(&mut out, opt, PrintState::new()).unwrap();
                if commas{
                    out.push(',');
                } else {
                    out.push(' ');
                }
                b.format(&mut out, opt, PrintState::new()).unwrap();
                if parens{
                    out.push(')');
                }
                Some(out)
            }
            _=>None}

    }),
    // sharp: symbol!("♯";Symmetric),
    metric: symbol!(METRIC_NAME;Symmetric,Real;print = |a, opt| {

        match opt.custom_print_mode {
             Some(("typst", 1)) =>{
                 if let AtomView::Fun(_)=a {
                     let body = r#"(a,b) = {
if a.at("lower",default:false) and b.at("lower",default:false){
$eta_(#to-eq(a) #to-eq(b))$
} else if a.at("lower",default:false) and b.at("upper",default:false){
$delta_(#to-eq(a))^#to-eq(b)$
}else if b.at("lower",default:false) and a.at("upper",default:false){
$delta_(#to-eq(b))^#to-eq(a)$
}else if a.at("upper",default:false) and b.at("upper",default:false){
$eta^(#to-eq(a)^#to-eq(b))$
} else{
$g(#to-eq(a),#to-eq(b))$
}
}"#;
                     return Some(body.into())
                 }
             }
             Some(("spenso",i))=>{
                 let SpensoPrintSettings{
                     parens,
                     commas,..
                 } = SpensoPrintSettings::from(i);
                let AtomView::Fun(f)=a else {
                    return None;
                };


                if f.get_nargs()==2 {
                    let mut argitem = f.iter();
                    let a = argitem.next().unwrap();
                    let b = argitem.next().unwrap();

                    let AtomView::Fun(mut f_a)=a else {
                     return None;
                    };
                    let AtomView::Fun(mut f_b)=b else {
                     return None;
                    };

                    let mut a_sym = f_a.get_symbol();
                    let mut b_sym = f_b.get_symbol();

                    let a_is_dind = a_sym == AIND_SYMBOLS.dind;
                    if a_is_dind {
                        if f_a.get_nargs()==1{
                            let new_a = f_a.iter().next().unwrap();
                            let AtomView::Fun(new_f_a)=new_a else {
                                return None;
                            };
                            a_sym = new_f_a.get_symbol();
                            f_a = new_f_a;
                        }else{
                            return None;
                        }
                    }
                    let b_is_dind = b_sym == AIND_SYMBOLS.dind;
                    if b_is_dind {
                        if f_b.get_nargs()==1{
                            let new_b = f_b.iter().next().unwrap();
                            let AtomView::Fun(new_f_b)=new_b else {
                                return None;
                            };
                            b_sym = new_f_b.get_symbol();
                            f_b = new_f_b;
                        }else{
                            return None;
                        }
                    }

                    if a_sym != b_sym {
                        return None;
                    }

                    match (a_is_dind,b_is_dind){
                        (true,true)=>{
                            let mut out = if opt.color_builtin_symbols {
                                nu_ansi_term::Color::Magenta.paint("g_").to_string()
                            } else {
                                "g_".to_string()
                            };

                            if parens{
                                out.push('(');
                            }
                            f_a.as_view().format(&mut out, opt, PrintState::new()).unwrap();

                            if commas{
                                out.push(',');
                            } else {
                                out.push(' ');
                            }
                            f_b.as_view().format(&mut out, opt, PrintState::new()).unwrap();
                            if parens{
                                out.push(')');
                            }
                            return Some(out)

                        }
                        (true,false)=>{
                            let mut out = if opt.color_builtin_symbols {
                                nu_ansi_term::Color::Magenta.paint("δ_").to_string()
                            } else {
                                "δ_".to_string()
                            };
                            if parens{
                                out.push('(');
                            }
                            f_a.as_view().format(&mut out, opt, PrintState::new()).unwrap();
                            if parens{
                                out.push(')');
                            }
                            if opt.color_builtin_symbols {
                                out.push_str( &nu_ansi_term::Color::Magenta.paint("^").to_string())
                            } else {
                                out.push('^');
                            };
                            if parens{
                                out.push('(');
                            }
                            f_b.as_view().format(&mut out, opt, PrintState::new()).unwrap();
                            if parens{
                                out.push(')');
                            }
                            return Some(out)
                        }
                        (false,true)=>{
                            let mut out = if opt.color_builtin_symbols {
                                nu_ansi_term::Color::Magenta.paint("δ_").to_string()
                            } else {
                                "δ_".to_string()
                            };
                            if parens{
                                out.push('(');
                            }
                            f_b.as_view().format(&mut out, opt, PrintState::new()).unwrap();
                            if parens{
                                out.push(')');
                            }
                            if opt.color_builtin_symbols {
                                out.push_str( &nu_ansi_term::Color::Magenta.paint("^").to_string())
                            } else {
                                out.push('^');
                            };
                            if parens{
                                out.push('(');
                            }
                            f_a.as_view().format(&mut out, opt, PrintState::new()).unwrap();
                            if parens{
                                out.push(')');
                            }
                            return Some(out)
                        }
                        (false,false)=>{
                            let mut out = if opt.color_builtin_symbols {
                                nu_ansi_term::Color::Magenta.paint("g^").to_string()
                            } else {
                                "g^".to_string()
                            };
                            if parens{
                                out.push('(');
                            }
                            f_a.as_view().format(&mut out, opt, PrintState::new()).unwrap();
                            if commas{
                                out.push(',');
                            } else {
                                out.push(' ');
                            }
                            f_b.as_view().format(&mut out, opt, PrintState::new()).unwrap();
                            if parens{
                                out.push(')');
                            }
                            return Some(out)
                        }
                    }
                }
            },
            _=>{}
        }
        None
    }),
});

impl IdentityName for Symbol {
    fn id() -> Self {
        ETS.metric
    }
}
impl<T: HasStructure<Structure = ExplicitKey<Aind>>, Aind: AbsInd> Default
    for TensorLibrary<T, Aind>
{
    fn default() -> Self {
        Self {
            explicit_dimension: AHashMap::new(),
            generic_dimension: AHashMap::new(),
        }
    }
}

impl<T: HasStructure<Structure = ExplicitKey<Aind>>, Aind: AbsInd> TensorLibrary<T, Aind> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn merge(&mut self, other: &mut Self) {
        self.explicit_dimension
            .extend(other.explicit_dimension.drain());
        self.generic_dimension
            .extend(other.generic_dimension.drain());
    }
}

impl<T: TensorLibraryData> TensorLibraryData for ConcreteOrParam<T> {
    fn one() -> Self {
        ConcreteOrParam::Concrete(T::one())
    }

    fn minus_one() -> Self {
        ConcreteOrParam::Concrete(T::minus_one())
    }

    fn zero() -> Self {
        ConcreteOrParam::Concrete(T::zero())
    }
}

impl<
    Aind: AbsInd,
    T: HasStructure<Structure = ExplicitKey<Aind>> + Clone,
    S: TensorStructure + HasName<Name: IntoSymbol, Args: IntoArgs>,
> Library<S> for TensorLibrary<T, Aind>
{
    type Key = ExplicitKey<Aind>;
    type Value = PermutedStructure<T>;
    // type Structure = ExplicitKey;

    fn get<'a>(&'a self, key: &Self::Key) -> Result<Cow<'a, Self::Value>, LibraryError<Self::Key>> {
        // println!("Trying:{}", key);
        if let Some(tensor) = self.explicit_dimension.get(key) {
            // println!("found explicit");
            Ok(Cow::Borrowed(tensor))
        } else if let Some(builder) = self.generic_dimension.get(&key.clone().into()) {
            let permutation = PermutedStructure {
                structure: builder(key.clone()),
                rep_permutation: Permutation::id(key.order()),
                index_permutation: Permutation::id(key.order()),
            };
            // println!("found generic");
            Ok(Cow::Owned(permutation))
        } else {
            Err(LibraryError::NotFound(key.clone()))
        }
    }

    fn key_for_structure(
        &self,
        structure: &PermutedStructure<S>,
    ) -> Result<Self::Key, LibraryError<Self::Key>>
    where
        S: TensorStructure,
    {
        let a = ExplicitKey::from_structure(structure);
        if let Some(key) = a {
            if <TensorLibrary<T, Aind> as Library<S>>::get(self, &key).is_ok() {
                Ok(key)
            } else {
                Err(LibraryError::InvalidKey)
            }
        } else {
            Err(LibraryError::InvalidKey)
        }
    }
}

impl<
    Aind: AbsInd,
    T: HasStructure<Structure = ExplicitKey<Aind>>
        + SetTensorData<SetData = <T as LibraryTensor>::Data>
        + Clone
        + LibraryTensor
        + PermuteTensor<Permuted = T>,
> TensorLibrary<T, Aind>
{
    pub fn contains_explicit_name(&self, name: Symbol) -> bool {
        self.explicit_dimension
            .keys()
            .any(|key| key.name().is_some_and(|key_name| key_name == name))
    }

    pub fn explicit_names(&self) -> Vec<Symbol> {
        self.explicit_dimension
            .keys()
            .filter_map(|key| key.name())
            .collect()
    }

    pub fn get_key_from_name(
        &self,
        name: Symbol,
    ) -> Result<ExplicitKey<Aind>, LibraryError<ExplicitKey<Aind>>> {
        let keys: Vec<_> = self
            .explicit_dimension
            .keys()
            .filter(|k| k.name().unwrap() == name)
            .collect();

        match keys.len() {
            0 => Err(LibraryError::InvalidKey),
            1 => Ok(keys[0].clone()),
            _ => Err(LibraryError::MultipleKeys(name.to_string())),
        }
    }
    pub fn metric_key(rep: LibraryRep) -> ExplicitKey<Aind> {
        ExplicitKey::from_iter([rep.new_rep(4), rep.new_rep(4)], ETS.metric, None).structure
    }

    pub fn generic_mink_metric(key: ExplicitKey<Aind>) -> T
    where
        T::SetData: TensorLibraryData,
    {
        let dim: usize = key.get_dim(0).unwrap().try_into().unwrap();
        let mut tensor = T::empty(key, T::SetData::zero());

        for i in 0..dim {
            if i > 0 {
                tensor.set(&[i, i], -T::SetData::one()).unwrap();
            } else {
                tensor.set(&[i, i], T::SetData::one()).unwrap();
            }
        }
        tensor
    }

    pub fn update_ids(&mut self)
    where
        T::SetData: TensorLibraryData,
    {
        initialize();
        for r in LibraryRep::all_dualizables() {
            self.insert_generic(Self::id(*r), Self::checked_identity);
        }

        for r in LibraryRep::all_self_duals() {
            let id_metric = GenericKey::new(ETS.metric, None, vec![*r, *r]);
            self.insert_generic(id_metric, Self::checked_identity);
        }

        for r in LibraryRep::all_inline_metrics() {
            self.insert_generic(
                GenericKey::new(ETS.flat, None, vec![*r, *r]),
                Self::checked_identity,
            );
            let id_metric = GenericKey::new(ETS.metric, None, vec![*r, *r]);
            self.insert_generic(id_metric, Self::diag_unimodular_metric);
        }
    }

    pub fn id(rep: LibraryRep) -> GenericKey {
        GenericKey::new(ETS.metric, None, vec![rep, rep.dual()])
    }

    pub fn checked_identity(key: ExplicitKey<Aind>) -> T
    where
        T::SetData: TensorLibraryData,
    {
        debug_assert!(key.order() == 2);
        debug_assert!(key.get_rep(0).map(|r| r.dual()) == key.get_rep(1));

        Self::identity(key)
    }

    pub fn identity(key: ExplicitKey<Aind>) -> T
    where
        T::SetData: TensorLibraryData,
    {
        let dim: usize = key.get_dim(0).unwrap().try_into().unwrap();
        let mut tensor = T::empty(key, T::SetData::zero());

        for i in 0..dim {
            tensor.set(&[i, i], T::SetData::one()).unwrap();
        }
        tensor
    }

    pub fn diag_unimodular_metric(key: ExplicitKey<Aind>) -> T
    where
        T::SetData: TensorLibraryData,
    {
        let dim: usize = key.get_dim(0).unwrap().try_into().unwrap();
        let rep = key.get_rep(0).unwrap();
        let mut tensor = T::empty(key, T::SetData::zero());

        for i in 0..dim {
            if rep.is_neg(i) {
                tensor.set(&[i, i], T::SetData::minus_one()).unwrap();
            } else {
                tensor.set(&[i, i], T::SetData::one()).unwrap();
            }
        }
        tensor
    }

    pub fn insert_explicit(&mut self, data: PermutedStructure<T>) {
        let key = data.structure.structure().clone();
        self.explicit_dimension.insert(key, data);
    }

    pub fn insert_explicit_dense(
        &mut self,
        key: PermutedStructure<ExplicitKey<Aind>>,
        data: Vec<T::Data>,
    ) -> Result<()> {
        let tensor = T::from_dense(key.structure.clone(), data)?;
        let perm_tensor = PermutedStructure {
            rep_permutation: key.rep_permutation,
            index_permutation: key.index_permutation,
            structure: tensor,
        };

        self.explicit_dimension
            .insert(key.structure, perm_tensor.permute_inds_wrapped());
        Ok(())
    }

    pub fn insert_explicit_sparse(
        &mut self,
        key: PermutedStructure<ExplicitKey<Aind>>,
        data: impl IntoIterator<Item = (Vec<ConcreteIndex>, T::Data)>,
        zero: T::Data,
    ) -> Result<()> {
        let tensor = T::from_sparse(key.structure.clone(), data, zero)?;

        let perm_tensor = PermutedStructure {
            rep_permutation: key.rep_permutation.clone(),
            index_permutation: key.index_permutation.clone(),
            structure: tensor,
        };

        self.explicit_dimension
            .insert(key.structure, perm_tensor.permute_inds_wrapped());
        Ok(())
    }

    pub fn insert_generic(&mut self, key: GenericKey, data: fn(ExplicitKey<Aind>) -> T) {
        self.generic_dimension.insert(key, data);
    }

    pub fn get(&self, key: &ExplicitKey<Aind>) -> Result<Cow<'_, T>>
    where
        T: Clone,
        LibraryError<ExplicitKey<Aind>>: Into<eyre::Error>,
    {
        // println!("Trying:{}", key);

        if let Some(tensor) = self.explicit_dimension.get(key) {
            Ok(Cow::Borrowed(&tensor.structure))
        } else if let Some(builder) = self.generic_dimension.get(&key.clone().into()) {
            Ok(Cow::Owned(builder(key.clone())))
        } else {
            Err(LibraryError::NotFound(key.clone()).into())
        }
    }
}

#[cfg(test)]
#[allow(clippy::type_complexity)]
mod test {
    use symbolica::{function, parse};

    use crate::{
        network::{
            ExecutionResult, Network, Sequential, SmallestDegree, TensorOrScalarOrKey,
            library::panicing::ErroringLibrary,
            parsing::{ParseSettings, ShadowedStructure},
            store::NetworkStore,
        },
        shadowing::Concretize,
        structure::{
            ToSymbolic,
            abstract_index::AbstractIndex,
            representation::{Euclidean, Minkowski},
        },
        tensors::data::SparseOrDense,
    };

    use super::*;

    #[test]
    fn add_to_lib() {
        let mut lib =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();
        let key = ExplicitKey::from_iter(
            [
                Euclidean {}.new_rep(4).cast(),
                Euclidean {}.new_rep(4).cast(),
                LibraryRep::from(Minkowski {}).new_rep(4),
            ],
            symbol!("gamma"),
            None,
        );

        println!("{}", key.structure);
        println!("{}", key.rep_permutation);

        let one = ConcreteOrParam::Concrete(RealOrComplex::Real(1.));
        lib.insert_explicit_sparse(
            (key).clone(),
            [(vec![0, 0, 1], one)],
            ConcreteOrParam::Concrete(RealOrComplex::Real(0.)),
        )
        .unwrap();

        lib.get(&key.structure).unwrap();
        let indexed = key.clone().reindex([0, 1, 2]).unwrap().structure;
        let expr = indexed.to_symbolic(None).unwrap();
        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(expr.as_view(), &lib, &ParseSettings::default())
        .unwrap();

        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );
    }
    #[test]
    fn libperm() {
        let mut lib =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();
        let key = ExplicitKey::from_iter(
            [
                Euclidean {}.new_rep(2).cast(),
                Euclidean {}.new_rep(2).cast(),
                LibraryRep::from(Minkowski {}).new_rep(2),
            ],
            symbol!("gamma"),
            None,
        );

        println!("{}", key.structure);
        println!("{}", key.rep_permutation);

        let tensor = MixedTensor::Param(
            ParamTensor::from_sparse(
                key.structure.clone(),
                [
                    (vec![0, 0, 0], parse!("a")),
                    (vec![0, 0, 1], parse!("b")),
                    (vec![0, 1, 0], parse!("c")),
                    (vec![0, 1, 1], parse!("d")),
                    (vec![1, 0, 0], parse!("e")),
                    (vec![1, 0, 1], parse!("f")),
                    (vec![1, 1, 0], parse!("g")),
                    (vec![1, 1, 1], parse!("h")),
                ],
                Atom::Zero,
            )
            .unwrap()
            .to_dense(),
        );

        lib.insert_explicit(PermutedStructure {
            structure: tensor,
            rep_permutation: key.rep_permutation.clone(),
            index_permutation: key.index_permutation.clone(),
        });

        lib.get(&key.structure).unwrap();
        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(
            parse!("gamma(euc(2,1),euc(2,2),mink(2,0))-gamma(mink(2,0),euc(2,2),euc(2,1))")
                .as_view(),
            &lib,
            &ParseSettings::default(),
        )
        .unwrap();

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();

        print!("One {}", net.result_tensor(&lib).unwrap());
        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(
            parse!("gamma(mink(2,0),euc(2,2),euc(2,1))").as_view(),
            &lib,
            &ParseSettings::default(),
        )
        .unwrap();

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();

        print!("Two{}", net.result_tensor(&lib).unwrap())
    }

    #[test]
    fn not_in_lib() {
        let lib =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();
        let key = ExplicitKey::from_iter(
            [
                Euclidean {}.new_rep(4).cast(),
                LibraryRep::from(Minkowski {}).new_rep(4),
                Euclidean {}.new_rep(4).cast(),
            ],
            symbol!("gamma"),
            None,
        );

        let indexed = key.reindex([0, 1, 2]).unwrap().structure;
        let expr = indexed.to_symbolic(None).unwrap();
        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(expr.as_view(), &lib, &ParseSettings::default())
        .unwrap();

        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );

        if let ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) =
            net.result().unwrap()
        {
            // println!("YaY:{a}");
            println!("{tensor}");
            assert_eq!(tensor, &indexed.to_shell().concretize(None));
        } else {
            panic!("Not Key")
        }
    }

    #[test]
    fn flat() {
        let mut lib =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();
        lib.update_ids();
        fn p(m: impl Into<AbstractIndex>) -> Atom {
            let m_atom: AbstractIndex = m.into();
            let m_atom: Atom = m_atom.into();
            let mink = Minkowski {}.new_rep(4);
            function!(symbol!("spenso::p"), mink.to_symbolic([m_atom]))
        }

        let mink = Minkowski {}.new_rep(4);

        let a = mink.flat(1, 2) * p(2);

        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(a.as_view(), &lib, &ParseSettings::default())
        .map_err(|a| a.to_string())
        .unwrap();

        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();

        if let Ok(ExecutionResult::Val(v)) = net.result_tensor(&lib) {
            println!("{}", v)
        }
    }

    #[test]
    fn dot() {
        let lib =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();

        let expr = parse!("p(1,mink(4,2))*q(2,mink(4,2))");
        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(expr.as_view(), &lib, &ParseSettings::default())
        .map_err(|a| a.to_string())
        .unwrap();

        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );

        if let Ok(ExecutionResult::Val(a)) = net.result_scalar() {
            if let ConcreteOrParam::Param(a) = a.as_ref() {
                let res = parse!(
                    "p(1,cind(0))*q(2,cind(0))-p(1,cind(1))*q(2,cind(1))-p(1,cind(2))*q(2,cind(2))-p(1,cind(3))*q(2,cind(3))"
                );
                assert_eq!(a, &res);
            } else {
                panic!("Not Key")
            }
        } else {
            panic!("Not Key")
        }
    }

    #[test]
    fn big_expr() {
        initialize();
        let mut lib =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();
        lib.update_ids();

        let expr = parse!(
            " -G^2*(-g(mink(4,5),mink(4,6))*Q(2,mink(4,7))+g(mink(4,5),mink(4,6))*Q(3,mink(4,7))+g(mink(4,5),mink(4,7))*Q(2,mink(4,6))+g(mink(4,5),mink(4,7))*Q(4,mink(4,6))-g(mink(4,6),mink(4,7))*Q(3,mink(4,5))-g(mink(4,6),mink(4,7))*Q(4,mink(4,5)))*g(mink(4,2),mink(4,5))*g(mink(4,3),mink(4,6))*g(euc(4,0),euc(4,5))*g(euc(4,1),euc(4,4))*g(mink(4,4),mink(4,7))*vbar(1,euc(4,1))*u(0,euc(4,0))*ϵbar(2,mink(4,2))*ϵbar(3,mink(4,3))*gamma(euc(4,5),euc(4,4),mink(4,4))"
        );
        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(expr.as_view(), &lib, &ParseSettings::default())
        .map_err(|a| a.to_string())
        .unwrap();

        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();

        if let Ok(ExecutionResult::Val(v)) = net.result_scalar() {
            println!("Hi{}", v)
        }
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );
    }

    #[test]
    fn small_expr() {
        initialize();
        let mut lib =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();
        lib.update_ids();

        let expr = parse!(
            "(-g(mink(4,5),mink(4,6))*Q(2,mink(4,7))+g(mink(4,5),mink(4,6))*Q(3,mink(4,7)))*g(mink(4,2),mink(4,5))*g(mink(4,3),mink(4,6))*g(mink(4,4),mink(4,7))*ϵbar(2,mink(4,2))*ϵbar(3,mink(4,3))"
        );
        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(expr.as_view(), &lib, &ParseSettings::default())
        .map_err(|a| a.to_string())
        .unwrap();

        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();

        if let Ok(ExecutionResult::Val(v)) = net.result_tensor(&lib) {
            println!("{}", v)
        }
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.name().map(|a| a.to_string()).unwrap_or("".to_owned()),
                |_| "".to_string()
            )
        );
    }

    #[test]
    fn one_times_x() {
        let mut a: Network<
            NetworkStore<MixedTensor<f64, ShadowedStructure<AbstractIndex>>, Atom>,
            DummyKey,
            DummyKey,
        > = Network::one() * Network::from_scalar(Atom::var(symbol!("x")));

        // a.merge_ops();
        a.execute::<Sequential, SmallestDegree, _, _, _>(
            &DummyLibrary::default(),
            &ErroringLibrary::new(),
        )
        .unwrap();

        let res = a.result_scalar();
        if let Ok(ExecutionResult::Val(v)) = res {
            println!("Hi{}", v)
        } else {
            // panic!("AAAAA{}", a.dot())
        }
    }

    #[test]
    fn transposition() {
        initialize();
        let mut lib =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();
        lib.update_ids();

        let key = ExplicitKey::from_iter(
            [Euclidean {}.new_rep(4), Euclidean {}.new_rep(4)],
            symbol!("A"),
            None,
        )
        .structure;

        let mut a: DataTensor<_, _> = DenseTensor::fill(key.clone(), Atom::num(1)).into();
        a.set(&[3, 0], parse!("a")).unwrap();
        let a =
            PermutedStructure::identity(MixedTensor::<f64, ExplicitKey<AbstractIndex>>::param(a));

        lib.insert_explicit(a);
        #[allow(non_snake_case)]
        fn A(i: impl Into<AbstractIndex>, j: impl Into<AbstractIndex>) -> Atom {
            let euc = Euclidean {}.new_rep(4);
            function!(symbol!("A"), euc.slot(i).to_atom(), euc.slot(j).to_atom())
        }

        let expr = A(0, 1) - A(1, 0);

        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(expr.as_view(), &lib, &ParseSettings::default())
        .map_err(|a| a.to_string())
        .unwrap();

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();

        if let Ok(ExecutionResult::Val(v)) = net.result_tensor(&lib) {
            println!("{}", v)
        }
    }
    #[test]
    fn trace_metric_kron() {
        initialize();
        let mut lib =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();
        lib.update_ids();

        let mink = Minkowski {}.new_rep(4);

        let expr = mink.id(1, 0) * mink.id(0, 1);

        let mut net = Network::<
            NetworkStore<
                MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
                ConcreteOrParam<RealOrComplex<f64>>,
            >,
            _,
            Symbol,
        >::try_from_view(expr.as_view(), &lib, &ParseSettings::default())
        .map_err(|a| a.to_string())
        .unwrap();

        let fnlib = ErroringLibrary::new();
        net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();

        if let Ok(ExecutionResult::Val(v)) = net.result_tensor(&lib) {
            println!("{}", v)
        }
    }
}
