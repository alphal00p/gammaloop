use crate::{
    shadowing::symbolica_utils::{IntoArgs, IntoSymbol},
    structure::{
        HasName, HasStructure, TensorShell, TensorStructure, ToSymbolic,
        concrete_index::FlatIndex,
        slot::{IsAbstractSlot, ParseableAind},
    },
    tensors::{
        data::{DataTensor, DenseTensor},
        parametric::{MixedTensor, ParamTensor, TensorCoefficient},
        // symbolic::SymbolicTensor,
    },
};
use eyre::Result;
use linnet::permutation::Permutation;
use symbolica::{atom::Atom, evaluate::FunctionMap};

pub mod symbolica_utils;

/// Trait that enables shadowing of a tensor
///
/// This creates a dense tensor of atoms, where the atoms are the expanded indices of the tensor, with the global name as the name of the labels.
pub trait Shadowable:
    HasStructure<
        Structure: TensorStructure + HasName<Name: IntoSymbol, Args: IntoArgs> + Clone + Sized,
    > + Sized
where
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    // type Const;
    fn expanded_shadow(&self) -> Result<DenseTensor<Atom, Self::Structure>> {
        self.shadow(Self::Structure::expanded_coef)
    }
    fn expanded_shadow_perm(
        &self,
        perm: &Permutation,
    ) -> Result<DenseTensor<Atom, Self::Structure>> {
        self.shadow(|a, id| Self::Structure::expanded_coef_perm(a, id, perm))
    }

    fn flat_shadow(&self) -> Result<DenseTensor<Atom, Self::Structure>> {
        self.shadow(Self::Structure::flat_coef)
    }

    fn shadow<T>(
        &self,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> T,
    ) -> Result<DenseTensor<Atom, Self::Structure>>
    where
        T: TensorCoefficient,
    {
        self.structure().clone().to_dense_labeled(index_to_atom)
    }
}

pub trait Concretize<T> {
    fn concretize(self, perm: Option<Permutation>) -> T;
}

impl<S: Shadowable> Concretize<DenseTensor<Atom, S::Structure>> for S
where
    <<S::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn concretize(self, perm: Option<Permutation>) -> DenseTensor<Atom, S::Structure> {
        // self.flat_s
        // todo!()
        if let Some(perm) = perm {
            self.expanded_shadow_perm(&perm).unwrap()
        } else {
            self.expanded_shadow().unwrap()
        }
    }
}

impl<S: Shadowable> Concretize<DataTensor<Atom, S::Structure>> for S
where
    <<S::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn concretize(self, perm: Option<Permutation>) -> DataTensor<Atom, S::Structure> {
        // self.flat_s
        // todo!()
        <S as Concretize<DenseTensor<Atom, S::Structure>>>::concretize(self, perm).into()
    }
}

impl<S: Shadowable> Concretize<ParamTensor<S::Structure>> for S
where
    <<S::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn concretize(self, perm: Option<Permutation>) -> ParamTensor<S::Structure> {
        // self.flat_s
        // todo!()
        ParamTensor::param(
            <S as Concretize<DataTensor<Atom, S::Structure>>>::concretize(self, perm),
        )
    }
}

impl<T: Clone, S: Shadowable> Concretize<MixedTensor<T, S::Structure>> for S
where
    <<S::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn concretize(self, perm: Option<Permutation>) -> MixedTensor<T, S::Structure> {
        // self.flat_s
        // todo!()
        MixedTensor::<T, S::Structure>::param(
            <S as Concretize<DataTensor<Atom, S::Structure>>>::concretize(self, perm),
        )
    }
}

// pub type ExplicitKey = IndexlessNamedStructure<Symbol, Vec<Atom>, Rep>;

// impl ExplicitKey {
//     pub fn from_structure<S: TensorStructure + HasName<Name: IntoSymbol, Args: IntoArgs>>(
//         structure: S,
//     ) -> Option<Self>
//     where
//         Rep: From<<S::Slot as IsAbstractSlot>::R>,
//     {
//         Some(IndexlessNamedStructure::from_iter(
//             structure.reps().into_iter().map(|r| r.cast()),
//             structure.name()?.ref_into_symbol(),
//             structure.args().map(|a| a.args()),
//         ))
//     }
// }

// #[derive(Clone, PartialEq, Eq, Debug, Hash)]
// pub struct GenericKey {
//     global_name: Symbol,
//     args: Option<Vec<Atom>>,
//     reps: Vec<Rep>,
// }

// impl GenericKey {
//     pub fn new(global_name: Symbol, args: Option<Vec<Atom>>, reps: Vec<Rep>) -> Self {
//         Self {
//             global_name,
//             args,
//             reps,
//         }
//     }
// }

// impl From<ExplicitKey> for GenericKey {
//     fn from(key: ExplicitKey) -> Self {
//         Self {
//             global_name: key.global_name.unwrap(),
//             reps: key.reps().into_iter().map(|r| r.rep).collect(),
//             args: key.additional_args,
//         }
//     }
// }

// #[allow(clippy::type_complexity)]
// pub struct ExplicitTensorMap<Data: Clone = f64> {
//     explicit_dimension: AHashMap<ExplicitKey, MixedTensor<Data, ExplicitKey>>,
//     generic_dimension: AHashMap<GenericKey, fn(ExplicitKey) -> MixedTensor<Data, ExplicitKey>>,
// }

// #[derive(Debug, Error)]
// pub enum ExplicitTensorError {
//     #[error("The key {0} is not present")]
//     KeyNotFound(ExplicitKey),
// }

// pub struct ExplicitTensorSymbols {
//     pub id: Symbol,
//     pub gamma: Symbol,
//     pub gamma5: Symbol,
//     pub proj_m: Symbol,
//     pub proj_p: Symbol,
//     pub sigma: Symbol,
//     pub metric: Symbol,
// }

// pub static ETS: LazyLock<ExplicitTensorSymbols> = LazyLock::new(|| ExplicitTensorSymbols {
//     id: symbol!(ExplicitTensorMap::<f64>::ID_NAME;Symmetric).unwrap(),
//     gamma: symbol!(ExplicitTensorMap::<f64>::GAMMA_NAME),
//     gamma5: symbol!(ExplicitTensorMap::<f64>::GAMMA5_NAME),
//     proj_m: symbol!(ExplicitTensorMap::<f64>::PROJ_M_NAME),
//     proj_p: symbol!(ExplicitTensorMap::<f64>::PROJ_P_NAME),
//     sigma: symbol!(ExplicitTensorMap::<f64>::SIGMA_NAME),
//     metric: symbol!(ExplicitTensorMap::<f64>::METRIC_NAME;Symmetric).unwrap(),
// });

// impl<Data: Clone> ExplicitTensorMap<Data> {
//     pub const ID_NAME: &'static str = "id";
//     pub const GAMMA_NAME: &'static str = "γ";
//     pub const GAMMA5_NAME: &'static str = "γ5";
//     pub const PROJ_M_NAME: &'static str = "ProjM";
//     pub const PROJ_P_NAME: &'static str = "ProjP";
//     pub const SIGMA_NAME: &'static str = "σ";
//     pub const METRIC_NAME: &'static str = "Metric";

//     pub fn new() -> Self
//     where
//         Data: num::One + num::Zero + Neg<Output = Data> + Default,
//     {
//         let mut new = Self {
//             explicit_dimension: AHashMap::new(),
//             generic_dimension: AHashMap::new(),
//         };
//         new.update_ids();

//         new.insert_generic_real(
//             Self::metric_key(ExtendibleReps::LORENTZ_UP).into(),
//             Self::generic_mink_metric,
//         );

//         new.insert_generic_real(
//             Self::metric_key(ExtendibleReps::MINKOWSKI).into(),
//             Self::generic_mink_metric,
//         );

//         new.insert_generic_real(
//             Self::metric_key(ExtendibleReps::LORENTZ_DOWN).into(),
//             Self::generic_mink_metric,
//         );

//         new.insert_generic_real(
//             Self::metric_key(ExtendibleReps::BISPINOR).into(),
//             Self::identity,
//         );

//         // println!("{:?}", Self::gamma_key());

//         new.insert_explicit_complex_sparse(Self::gamma_key(), Self::gamma_data_weyl_transposed());
//         new.insert_explicit_complex_sparse(Self::gamma_five_key(), Self::gamma_five_data_weyl());
//         new.insert_explicit_real_sparse(Self::proj_m_key(), Self::proj_m_data_weyl());
//         new.insert_explicit_real_sparse(Self::proj_p_key(), Self::proj_p_data_weyl());
//         new.insert_explicit_complex_sparse(Self::sigma_key(), Self::sigma_data_weyl());
//         new
//     }

//     // Gamma(1,2,3) Dirac matrix (γ^μ1)_s2_s3
//     fn gamma_key() -> ExplicitKey {
//         ExplicitKey::from_iter(
//             [
//                 ExtendibleReps::MINKOWSKI.new_rep(4),
//                 ExtendibleReps::BISPINOR.new_rep(4),
//                 ExtendibleReps::BISPINOR.new_rep(4),
//             ],
//             symbol!(Self::GAMMA_NAME),
//             None,
//         )
//     }

//     fn metric_key(rep: Rep) -> ExplicitKey {
//         ExplicitKey::from_iter(
//             [rep.new_rep(4), rep.new_rep(4)],
//             symbol!(Self::METRIC_NAME),
//             None,
//         )
//     }

//     fn generic_mink_metric(key: ExplicitKey) -> MixedTensor<Data, ExplicitKey>
//     where
//         Data: num::One + Neg<Output = Data> + Default,
//     {
//         let dim: usize = key.get_dim(0).unwrap().try_into().unwrap();
//         let mut tensor = SparseTensor::empty(key);

//         for i in 0..dim {
//             if i > 0 {
//                 tensor.set(&[i, i], -Data::one()).unwrap();
//             } else {
//                 tensor.set(&[i, i], Data::one()).unwrap();
//             }
//         }
//         tensor.into()
//     }

//     fn gamma_data_weyl_transposed() -> Vec<(Vec<ConcreteIndex>, Complex<Data>)>
//     where
//         Data: num::Zero + num::One + Neg<Output = Data> + Clone,
//     {
//         let c1 = Complex::<Data>::new(Data::one(), Data::zero());
//         let cn1 = Complex::<Data>::new(-Data::one(), Data::zero());
//         let ci = Complex::<Data>::new(Data::zero(), Data::one());
//         let cni = Complex::<Data>::new(Data::zero(), -Data::one());
//         vec![
//             (vec![0, 0, 2], c1.clone()),
//             (vec![0, 1, 3], c1.clone()),
//             (vec![0, 2, 0], c1.clone()),
//             (vec![0, 3, 1], c1.clone()),
//             (vec![1, 0, 3], c1.clone()),
//             (vec![1, 1, 2], c1.clone()),
//             (vec![1, 2, 1], cn1.clone()),
//             (vec![1, 3, 0], cn1.clone()),
//             (vec![2, 0, 3], cni.clone()),
//             (vec![2, 1, 2], ci.clone()),
//             (vec![2, 2, 1], ci.clone()),
//             (vec![2, 3, 0], cni.clone()),
//             (vec![3, 0, 2], c1.clone()),
//             (vec![3, 1, 3], cn1.clone()),
//             (vec![3, 2, 0], cn1.clone()),
//             (vec![3, 3, 1], c1.clone()),
//         ]
//     }

//     fn _gamma_data_weyl() -> Vec<(Vec<ConcreteIndex>, Complex<Data>)>
//     where
//         Data: num::Zero + num::One + Neg<Output = Data> + Clone,
//     {
//         let c1 = Complex::<Data>::new(Data::one(), Data::zero());
//         let cn1 = Complex::<Data>::new(-Data::one(), Data::zero());
//         let ci = Complex::<Data>::new(Data::zero(), Data::one());
//         let cni = Complex::<Data>::new(Data::zero(), -Data::one());
//         vec![
//             (vec![0, 2, 0], c1.clone()),
//             (vec![0, 3, 1], c1.clone()),
//             (vec![0, 0, 2], c1.clone()),
//             (vec![0, 1, 3], c1.clone()),
//             (vec![1, 3, 0], c1.clone()),
//             (vec![1, 2, 1], c1.clone()),
//             (vec![1, 1, 2], cn1.clone()),
//             (vec![1, 0, 3], cn1.clone()),
//             (vec![2, 3, 0], cni.clone()),
//             (vec![2, 2, 1], ci.clone()),
//             (vec![2, 1, 2], ci.clone()),
//             (vec![2, 0, 3], cni.clone()),
//             (vec![3, 2, 0], c1.clone()),
//             (vec![3, 3, 1], cn1.clone()),
//             (vec![3, 0, 2], cn1.clone()),
//             (vec![3, 1, 3], c1.clone()),
//         ]
//     }

//     fn gamma_five_key() -> ExplicitKey {
//         ExplicitKey::from_iter(
//             [
//                 ExtendibleReps::BISPINOR.new_rep(4),
//                 ExtendibleReps::BISPINOR.new_rep(4),
//             ],
//             symbol!(Self::GAMMA5_NAME),
//             None,
//         )
//     }

//     fn gamma_five_data_weyl() -> Vec<(Vec<ConcreteIndex>, Complex<Data>)>
//     where
//         Data: num::Zero + num::One + Neg<Output = Data> + Clone,
//     {
//         let c1 = Complex::<Data>::new(Data::one(), Data::zero());

//         vec![
//             (vec![0, 0], -c1.clone()),
//             (vec![1, 1], -c1.clone()),
//             (vec![2, 2], c1.clone()),
//             (vec![3, 3], c1.clone()),
//         ]
//     }

//     fn proj_m_key() -> ExplicitKey {
//         ExplicitKey::from_iter(
//             [
//                 ExtendibleReps::BISPINOR.new_rep(4),
//                 ExtendibleReps::BISPINOR.new_rep(4),
//             ],
//             symbol!(Self::PROJ_M_NAME),
//             None,
//         )
//     }

//     fn proj_p_key() -> ExplicitKey {
//         ExplicitKey::from_iter(
//             [
//                 ExtendibleReps::BISPINOR.new_rep(4),
//                 ExtendibleReps::BISPINOR.new_rep(4),
//             ],
//             symbol!(Self::PROJ_P_NAME),
//             None,
//         )
//     }

//     fn proj_m_data_weyl() -> Vec<(Vec<ConcreteIndex>, Data)>
//     where
//         Data: num::One,
//     {
//         vec![(vec![0, 0], Data::one()), (vec![1, 1], Data::one())]
//     }

//     fn proj_p_data_weyl() -> Vec<(Vec<ConcreteIndex>, Data)>
//     where
//         Data: num::One,
//     {
//         vec![(vec![2, 2], Data::one()), (vec![3, 3], Data::one())]
//     }

//     fn sigma_key() -> ExplicitKey {
//         ExplicitKey::from_iter(
//             [
//                 ExtendibleReps::MINKOWSKI.new_rep(4),
//                 ExtendibleReps::MINKOWSKI.new_rep(4),
//                 ExtendibleReps::BISPINOR.new_rep(4),
//                 ExtendibleReps::BISPINOR.new_rep(4),
//             ],
//             symbol!(Self::SIGMA_NAME),
//             None,
//         )
//     }

//     fn sigma_data_weyl() -> Vec<(Vec<ConcreteIndex>, Complex<Data>)>
//     where
//         Data: num::Zero + num::One + Neg<Output = Data> + Clone,
//     {
//         let c1 = Complex::<Data>::new(Data::one(), Data::zero());
//         let cn1 = Complex::<Data>::new(-Data::one(), Data::zero());
//         let ci = Complex::<Data>::new(Data::zero(), Data::one());
//         let cni = Complex::<Data>::new(Data::zero(), -Data::one());

//         vec![
//             (vec![0, 2, 0, 1], c1.clone()),
//             (vec![0, 2, 3, 0], c1.clone()),
//             (vec![0, 3, 1, 2], c1.clone()),
//             (vec![1, 0, 2, 2], c1.clone()),
//             (vec![1, 1, 1, 2], c1.clone()),
//             (vec![1, 3, 0, 2], c1.clone()),
//             (vec![2, 2, 1, 0], c1.clone()),
//             (vec![2, 2, 2, 1], c1.clone()),
//             (vec![2, 3, 3, 2], c1.clone()),
//             (vec![3, 0, 0, 2], c1.clone()),
//             (vec![3, 3, 2, 2], c1.clone()),
//             (vec![3, 1, 3, 2], c1.clone()),
//             (vec![0, 1, 3, 0], ci.clone()),
//             (vec![0, 3, 1, 1], ci.clone()),
//             (vec![0, 3, 2, 0], ci.clone()),
//             (vec![1, 0, 3, 3], ci.clone()),
//             (vec![1, 1, 0, 3], ci.clone()),
//             (vec![1, 1, 2, 0], ci.clone()),
//             (vec![2, 1, 1, 0], ci.clone()),
//             (vec![2, 3, 0, 0], ci.clone()),
//             (vec![2, 3, 3, 1], ci.clone()),
//             (vec![3, 0, 1, 3], ci.clone()),
//             (vec![3, 1, 0, 0], ci.clone()),
//             (vec![3, 1, 2, 3], ci.clone()),
//             (vec![0, 0, 3, 2], cn1.clone()),
//             (vec![0, 1, 0, 2], cn1.clone()),
//             (vec![0, 2, 1, 3], cn1.clone()),
//             (vec![1, 2, 0, 3], cn1.clone()),
//             (vec![1, 2, 1, 1], cn1.clone()),
//             (vec![1, 2, 2, 0], cn1.clone()),
//             (vec![2, 0, 1, 2], cn1.clone()),
//             (vec![2, 1, 2, 2], cn1.clone()),
//             (vec![2, 2, 3, 3], cn1.clone()),
//             (vec![3, 2, 0, 0], cn1.clone()),
//             (vec![3, 2, 2, 3], cn1.clone()),
//             (vec![3, 2, 3, 1], cn1.clone()),
//             (vec![0, 0, 2, 3], cni.clone()),
//             (vec![0, 0, 3, 1], cni.clone()),
//             (vec![0, 1, 1, 3], cni.clone()),
//             (vec![1, 0, 2, 1], cni.clone()),
//             (vec![1, 3, 0, 1], cni.clone()),
//             (vec![1, 3, 3, 0], cni.clone()),
//             (vec![2, 0, 0, 3], cni.clone()),
//             (vec![2, 0, 1, 1], cni.clone()),
//             (vec![2, 1, 3, 3], cni.clone()),
//             (vec![3, 0, 0, 1], cni.clone()),
//             (vec![3, 3, 1, 0], cni.clone()),
//             (vec![3, 3, 2, 1], cni.clone()),
//         ]
//     }

//     pub fn update_ids(&mut self)
//     where
//         Data: num::One + Default,
//     {
//         for rep in REPS.read().unwrap().reps() {
//             self.insert_generic_real(Self::id(*rep), Self::checked_identity);
//             let id_metric =
//                 GenericKey::new(symbol!(Self::METRIC_NAME), None, vec![*rep, rep.dual()]);
//             self.insert_generic_real(id_metric, Self::checked_identity);
//             if rep.dual() != *rep {
//                 self.insert_generic_real(Self::id(rep.dual()), Self::checked_identity);
//                 let id_metric =
//                     GenericKey::new(symbol!(Self::METRIC_NAME), None, vec![rep.dual(), *rep]);
//                 self.insert_generic_real(id_metric, Self::checked_identity);
//             }
//         }
//     }

//     pub fn id(rep: Rep) -> GenericKey {
//         GenericKey::new(symbol!(Self::ID_NAME), None, vec![rep, rep.dual()])
//     }

//     pub fn checked_identity(key: ExplicitKey) -> MixedTensor<Data, ExplicitKey>
//     where
//         Data: num::One + Default,
//     {
//         assert!(key.order() == 2);
//         assert!(key.get_rep(0).map(|r| r.dual()) == key.get_rep(1));

//         Self::identity(key)
//     }

//     pub fn identity(key: ExplicitKey) -> MixedTensor<Data, ExplicitKey>
//     where
//         Data: num::One + Default,
//     {
//         let dim: usize = key.get_dim(0).unwrap().try_into().unwrap();
//         let mut tensor = SparseTensor::empty(key);

//         for i in 0..dim {
//             tensor.set(&[i, i], Data::one()).unwrap();
//         }
//         tensor.to_dense().into()
//     }

//     pub fn insert_explicit(&mut self, data: MixedTensor<Data, ExplicitKey>) {
//         let key = data.structure().clone();
//         self.explicit_dimension.insert(key, data);
//     }

//     pub fn insert_explicit_real(&mut self, key: ExplicitKey, data: Vec<Data>) {
//         let tensor: MixedTensor<Data, ExplicitKey> =
//             DenseTensor::from_data(data, key.clone()).unwrap().into();
//         self.explicit_dimension.insert(key, tensor);
//     }

//     pub fn insert_explicit_complex(&mut self, key: ExplicitKey, data: Vec<Complex<Data>>) {
//         let tensor: MixedTensor<Data, ExplicitKey> =
//             DenseTensor::from_data(data, key.clone()).unwrap().into();
//         self.explicit_dimension.insert(key, tensor);
//     }

//     pub fn insert_explicit_real_sparse(
//         &mut self,
//         key: ExplicitKey,
//         data: impl IntoIterator<Item = (Vec<ConcreteIndex>, Data)>,
//     ) {
//         let tensor: MixedTensor<Data, ExplicitKey> =
//             SparseTensor::from_data(data, key.clone()).unwrap().into();
//         self.explicit_dimension.insert(key, tensor);
//     }

//     pub fn insert_explicit_complex_sparse(
//         &mut self,
//         key: ExplicitKey,
//         data: impl IntoIterator<Item = (Vec<ConcreteIndex>, Complex<Data>)>,
//     ) {
//         let tensor: MixedTensor<Data, ExplicitKey> =
//             SparseTensor::from_data(data, key.clone()).unwrap().into();
//         self.explicit_dimension.insert(key, tensor);
//     }

//     pub fn insert_generic_real(
//         &mut self,
//         key: GenericKey,
//         data: fn(ExplicitKey) -> MixedTensor<Data, ExplicitKey>,
//     ) {
//         self.generic_dimension.insert(key, data);
//     }

//     pub fn get(&self, key: &ExplicitKey) -> Result<MixedTensor<Data, ExplicitKey>> {
//         if let Some(tensor) = self.explicit_dimension.get(key) {
//             Ok(tensor.clone())
//         } else if let Some(builder) = self.generic_dimension.get(&key.clone().into()) {
//             Ok(builder(key.clone()))
//         } else {
//             Err(ExplicitTensorError::KeyNotFound(key.clone()).into())
//         }
//     }
// }

// impl<Data: Clone + num::One + num::Zero + Neg<Output = Data> + Default> Default
//     for ExplicitTensorMap<Data>
// {
//     fn default() -> Self {
//         Self::new()
//     }
// }

pub trait ShadowMapping<Const>: Shadowable
where
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn expanded_shadow_with_map(
        &self,
        fn_map: &mut FunctionMap<Const>,
    ) -> Result<ParamTensor<Self::Structure>> {
        self.shadow_with_map(fn_map, Self::Structure::expanded_coef)
    }

    fn shadow_with_map<T, F>(
        &self,
        fn_map: &mut FunctionMap<Const>,
        index_to_atom: F,
    ) -> Result<ParamTensor<Self::Structure>>
    where
        T: TensorCoefficient,
        F: Fn(&Self::Structure, FlatIndex) -> T + Clone,
    {
        // Some(ParamTensor::Param(self.shadow(index_to_atom)?.into()))
        self.append_map(fn_map, index_to_atom.clone());
        self.shadow(index_to_atom)
            .map(|x| ParamTensor::param(x.into()))
    }

    fn append_map<T>(
        &self,
        fn_map: &mut FunctionMap<Const>,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> T,
    ) where
        T: TensorCoefficient;

    fn flat_append_map(&self, fn_map: &mut FunctionMap<Const>) {
        self.append_map(fn_map, Self::Structure::flat_coef)
    }

    fn expanded_append_map(&self, fn_map: &mut FunctionMap<Const>) {
        self.append_map(fn_map, Self::Structure::expanded_coef)
    }

    fn flat_shadow_with_map(
        &self,
        fn_map: &mut FunctionMap<Const>,
    ) -> Result<ParamTensor<Self::Structure>> {
        self.shadow_with_map(fn_map, Self::Structure::flat_coef)
    }
}

impl<S: TensorStructure + HasName + Clone> Shadowable for TensorShell<S>
where
    S::Name: IntoSymbol + Clone,
    S::Args: IntoArgs,
    <<S as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
}

impl<S: TensorStructure + HasName + Clone, Const> ShadowMapping<Const> for TensorShell<S>
where
    S::Name: IntoSymbol + Clone,
    S::Args: IntoArgs,
    <<S as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn append_map<T>(
        &self,
        _fn_map: &mut FunctionMap<Const>,
        _index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> T,
    ) where
        T: TensorCoefficient,
    {
    }
}

#[cfg(test)]
pub mod test {
    use std::sync::RwLock;

    use once_cell::sync::Lazy;
    use symbolica::atom::Symbol;

    use crate::network::library::{symbolic::ETS, symbolic::ExplicitKey, symbolic::TensorLibrary};

    #[allow(clippy::type_complexity)]
    pub static EXPLICIT_TENSOR_MAP: Lazy<
        RwLock<TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>>,
    > = Lazy::new(|| {
        let mut lib = TensorLibrary::new();
        lib.update_ids();
        RwLock::new(lib)
    });

    use crate::structure::abstract_index::AbstractIndex;
    use crate::structure::{OrderedStructure, PermutedStructure};
    use crate::{
        contraction::Contract,
        structure::{
            HasStructure, IndexlessNamedStructure, TensorStructure,
            representation::{LibraryRep, RepName},
        },
        tensors::parametric::MixedTensor,
    };

    #[test]
    fn test_identity() {
        // EXPLICIT_TENSOR_MAP.write().unwrap().update_ids();
        let mut tensor_library: TensorLibrary<
            MixedTensor<f64, ExplicitKey<AbstractIndex>>,
            AbstractIndex,
        > = TensorLibrary::new();
        tensor_library.update_ids();

        for rep in LibraryRep::all_representations() {
            let structure = [rep.new_rep(4), rep.new_rep(4).dual()];

            let idstructure: PermutedStructure<IndexlessNamedStructure<Symbol, (), LibraryRep>> =
                IndexlessNamedStructure::from_iter(structure, ETS.metric, None);

            let idkey = ExplicitKey::from_structure(&idstructure).unwrap();

            let id = tensor_library.get(&idkey).unwrap().into_owned();

            let trace_structure: OrderedStructure = PermutedStructure::from_iter([
                rep.new_rep(4).slot(3),
                rep.new_rep(4).dual().slot(4),
            ])
            .structure;
            let id1 = id.map_structure(|_| trace_structure.clone());
            let id2 = id1
                .clone()
                .map_structure(|_| trace_structure.clone().dual());

            // println!("{}", rep);
            assert_eq!(
                4.,
                id1.contract(&id2)
                    .unwrap()
                    .scalar()
                    .unwrap()
                    .try_into_concrete()
                    .unwrap()
                    .try_into_real()
                    .unwrap(),
                "trace of 4-dim identity should be 4 for rep {}",
                rep
            );
        }
    }

    // #[test]
    // fn pslash() {
    //     let expr = parse!("p(1,mink(4,mu))*γ(mink(4,mu),bis(4,i),bis(4,j))").unwrap();
    //     let mut network: TensorNetwork<MixedTensor<f64, ShadowedStructure>, Atom> =
    //         TensorNetwork::<MixedTensor<f64, ShadowedStructure>, Atom>::try_from_view(
    //             expr.as_view(),
    //             &EXPLICIT_TENSOR_MAP.read().unwrap(),
    //         )
    //         .unwrap();
    //     network.contract().unwrap();
    //     let result = network
    //         .result()
    //         .unwrap()
    //         .0
    //         .try_into_parametric()
    //         .unwrap()
    //         .tensor
    //         .map_data(|a| a.to_string())
    //         .map_structure(VecStructure::from);

    //     insta::assert_ron_snapshot!(result);
    // }
}
