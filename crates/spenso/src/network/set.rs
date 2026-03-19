use serde::{Deserialize, Serialize};

use crate::{
    algebra::algebraic_traits::{One, Zero},
    structure::{
        PermutedStructure,
        permuted::PermuteTensor,
        slot::{AbsInd, IsAbstractSlot},
    },
};
use std::{
    borrow::Cow,
    fmt::{Debug, Display},
};

#[cfg(feature = "shadowing")]
use eyre::eyre;

#[cfg(feature = "shadowing")]
use symbolica::{
    atom::{Atom, AtomView},
    domains::{
        InternalOrdering,
        float::{Complex as SymComplex, Real, SingleFloat},
        rational::Rational,
    },
    evaluate::{
        CompileOptions, CompiledCode, CompiledNumber, EvalTree, ExportNumber, ExportSettings,
        ExportedCode, ExpressionEvaluator, FunctionMap, OptimizationSettings,
    },
};

#[cfg(feature = "shadowing")]
use crate::{
    algebra::complex::{Complex, symbolica_traits::CompiledComplexEvaluatorSpenso},
    tensors::data::{DataIterator, DenseTensor, SetTensorData, SparseTensor},
    tensors::parametric::ParamTensor,
};

#[cfg(feature = "shadowing")]
use super::store::TensorScalarStoreMapping;

use crate::{
    structure::{CastStructure, HasStructure, ScalarTensor, TensorStructure},
    tensors::data::DataTensor,
};

use super::{
    ExecutionResult,
    library::LibraryTensor,
    store::{NetworkStore, TensorScalarStore},
};
use super::{Library, Network, TensorNetworkError, TensorOrScalarOrKey};

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    // bincode_trait_derive::BorrowDecodeFromDecode,
)]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
pub struct TensorNetworkSet<S, K, FK, Aind> {
    pub networks: Vec<Network<S, K, FK, Aind>>,
}

impl<S, K, FK, Aind> FromIterator<Network<S, K, FK, Aind>> for TensorNetworkSet<S, K, FK, Aind> {
    fn from_iter<T: IntoIterator<Item = Network<S, K, FK, Aind>>>(iter: T) -> Self {
        Self {
            networks: Vec::from_iter(iter),
        }
    }
}

impl<S, K, FK, Aind> TensorNetworkSet<S, K, FK, Aind> {
    pub fn new() -> Self {
        TensorNetworkSet { networks: vec![] }
    }

    pub fn push(&mut self, network: Network<S, K, FK, Aind>) {
        // self.scalars.push(network.scalar);
        self.networks.push(network);
    }
}

impl<S, K, FK, Aind> Default for TensorNetworkSet<S, K, FK, Aind> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(feature = "shadowing")]
pub type EvalTreeTensorNetworkSet<T, S, K, FK, Aind, Str> =
    SharedTensorNetworkSet<EvalTree<T>, S, K, FK, Aind, Str>;

#[cfg(feature = "shadowing")]
pub type EvalTensorNetworkSet<T, S, K, FK, Aind, Str> =
    SharedTensorNetworkSet<ExpressionEvaluator<T>, S, K, FK, Aind, Str>;

#[cfg(feature = "shadowing")]
pub type CompiledTensorNetworkSet<S, K, FK, Aind, Str> =
    SharedTensorNetworkSet<CompiledComplexEvaluatorSpenso, S, K, FK, Aind, Str>;

#[derive(Debug, Clone)]
pub struct SharedTensorNetworkSet<
    D,
    S: TensorStructure,
    K,
    FK,
    Aind,
    Str: TensorScalarStore<Tensor = DataTensor<usize, S>, Scalar = usize> = NetworkStore<
        DataTensor<usize, S>,
        usize,
    >,
> {
    pub networks: Vec<Network<Str, K, FK, Aind>>,
    pub shared_data: D,
    pub len: usize,
}

impl<
    T: TensorStructure,
    S,
    K: Display + Debug,
    FK: Display + Debug + Clone,
    Str: TensorScalarStore<Tensor = T, Scalar = S>,
    Aind: AbsInd,
> TensorNetworkSet<Str, K, FK, Aind>
{
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn result(
        &self,
    ) -> Result<
        Vec<ExecutionResult<TensorOrScalarOrKey<&T, &S, &PermutedStructure<K>, Aind>>>,
        TensorNetworkError<K, FK>,
    >
    where
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        self.networks.iter().map(|n| n.result()).collect()
    }

    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn result_tensor<'a, LT, L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>>(
        &'a self,
        lib: &L,
    ) -> Result<Vec<ExecutionResult<Cow<'a, T>>>, TensorNetworkError<K, FK>>
    where
        S: 'a,
        T: Clone + ScalarTensor + HasStructure,
        T::Scalar: One + Zero,
        for<'b> &'b S: Into<T::Scalar>,
        LT: TensorStructure<Indexed = T> + Clone + LibraryTensor<WithIndices = T>,
        T: PermuteTensor<Permuted = T>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        self.networks.iter().map(|n| n.result_tensor(lib)).collect()
    }

    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn result_scalar<'a>(
        &'a self,
    ) -> Result<Vec<ExecutionResult<Cow<'a, S>>>, TensorNetworkError<K, FK>>
    where
        T: Clone + ScalarTensor + 'a,
        T::Scalar: Into<S>,
        S: One + Zero + Clone,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        self.networks.iter().map(|n| n.result_scalar()).collect()
    }

    pub fn cast<U>(self) -> TensorNetworkSet<Str::Store<U, S>, K, FK, Aind>
    where
        K: Clone,
        T: CastStructure<U> + HasStructure,
        T::Structure: TensorStructure,
        U: HasStructure,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        U::Structure: From<T::Structure> + TensorStructure<Slot = T::Slot>,
    {
        self.networks.into_iter().map(|n| n.cast()).collect()
    }
}

#[cfg(feature = "shadowing")]
impl<
    Store: TensorScalarStore<Tensor = ParamTensor<S>, Scalar = Atom>,
    S: TensorStructure + Clone,
    K: Clone,
    FK: Clone,
    Aind: AbsInd,
> TensorNetworkSet<Store, K, FK, Aind>
{
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn eval_tree(
        self,
        fn_map: &FunctionMap,
        params: &[Atom],
    ) -> eyre::Result<
        EvalTreeTensorNetworkSet<
            SymComplex<Rational>,
            S,
            K,
            FK,
            Aind,
            Store::Store<DataTensor<usize, S>, usize>,
        >,
    >
    where
        Store::Store<DataTensor<usize, S>, usize>:
            TensorScalarStore<Tensor = DataTensor<usize, S>, Scalar = usize>,
        S: TensorStructure,
    {
        let mut networks = vec![];

        let mut atoms = vec![];
        let mut scalar_id = 0;
        let mut tensor_id = 0;

        for net in &self.networks {
            // let refnet = net.map_ref(|a|, |t|t);
            for s in net.iter_scalars() {
                atoms.push(s);
                tensor_id += 1;
            }

            let tensors = net.map_ref(
                |_| {
                    let oldid = scalar_id;
                    scalar_id += 1;
                    oldid
                },
                |t| {
                    let structure = t.structure().clone();
                    match &t.tensor {
                        DataTensor::Dense(d) => {
                            let oldid = tensor_id;
                            tensor_id += d.size().unwrap();
                            for (_, a) in d.flat_iter() {
                                atoms.push(a);
                            }
                            DataTensor::Dense(
                                DenseTensor::from_data(Vec::from_iter(oldid..tensor_id), structure)
                                    .expect("Failed to create DenseTensor"),
                            )
                        }
                        DataTensor::Sparse(s) => {
                            let mut t = SparseTensor::empty(structure, 0);
                            for (i, a) in s.flat_iter() {
                                t.set_flat(i, tensor_id)
                                    .expect("Failed to set value in SparseTensor");
                                atoms.push(a);
                                tensor_id += 1;
                            }
                            DataTensor::Sparse(t)
                        }
                    }
                },
            );

            networks.push(tensors);
        }

        Ok(EvalTreeTensorNetworkSet {
            networks,
            shared_data: AtomView::to_eval_tree_multiple(&atoms, fn_map, params)
                .map_err(|s| eyre!(s))?,
            len: atoms.len(),
        })
    }
}

#[cfg(feature = "shadowing")]
impl<
    S: Clone + TensorStructure,
    K,
    FK,
    Aind: AbsInd,
    Store: TensorScalarStore<Tensor = DataTensor<usize, S>, Scalar = usize>,
> EvalTreeTensorNetworkSet<SymComplex<Rational>, S, K, FK, Aind, Store>
{
    pub fn horner_scheme(&mut self) {
        self.shared_data.horner_scheme();
    }
}
#[cfg(feature = "shadowing")]
impl<
    S: TensorStructure + Clone,
    K: Clone,
    FK: Clone,
    Aind: AbsInd,
    Store: TensorScalarStore<Tensor = DataTensor<usize, S>, Scalar = usize> + Clone,
> EvalTreeTensorNetworkSet<SymComplex<Rational>, S, K, FK, Aind, Store>
{
    pub fn linearize(
        mut self,
        settings: &OptimizationSettings,
    ) -> EvalTensorNetworkSet<SymComplex<Rational>, S, K, FK, Aind, Store> {
        EvalTensorNetworkSet {
            networks: self.networks,
            shared_data: self.shared_data.linearize(settings),
            len: self.len,
        }
    }
}

#[cfg(feature = "shadowing")]
impl<
    T,
    S: TensorStructure + Clone,
    K: Clone,
    FK: Clone,
    Aind: AbsInd,
    Store: TensorScalarStore<Tensor = DataTensor<usize, S>, Scalar = usize> + Clone,
> EvalTreeTensorNetworkSet<T, S, K, FK, Aind, Store>
{
    pub fn map_coeff<T2, F: Fn(&T) -> T2>(
        &self,
        f: &F,
    ) -> EvalTreeTensorNetworkSet<T2, S, K, FK, Aind, Store>
    where
        T: Clone + PartialEq,
    {
        EvalTreeTensorNetworkSet {
            networks: self.networks.clone(),
            shared_data: self.shared_data.map_coeff(f),
            len: self.len,
        }
        // self.map_data_ref(|x| x.map_coeff(f))
    }

    pub fn common_subexpression_elimination(&mut self)
    where
        T: std::fmt::Debug + std::hash::Hash + Eq + InternalOrdering + Clone + Default,
    {
        self.shared_data.common_subexpression_elimination()
    }

    #[allow(clippy::type_complexity)]
    pub fn evaluate(
        &mut self,
        params: &[T],
    ) -> TensorNetworkSet<Store::Store<DataTensor<T, S>, T>, K, FK, Aind>
    where
        T: Real + SingleFloat,
        S: TensorStructure + Clone,
    {
        let zero = params[0].zero();
        let mut data = vec![zero; self.len];

        let mut networks = vec![];

        self.shared_data.evaluate(params, &mut data);

        for net in self.networks.iter() {
            let data_net = net.map_ref(
                |s| data[*s].clone(),
                |p| {
                    let structure = p.structure().clone();
                    match &p {
                        DataTensor::Dense(d) => {
                            let mut t_data = vec![];
                            for (_, &a) in d.flat_iter() {
                                t_data.push(data[a].clone());
                            }
                            DataTensor::Dense(DenseTensor::from_data(t_data, structure).unwrap())
                        }
                        DataTensor::Sparse(s) => {
                            let mut t = SparseTensor::empty(structure, T::new_zero());
                            for (i, &a) in s.flat_iter() {
                                t.set_flat(i, data[a].clone()).unwrap();
                            }
                            DataTensor::Sparse(t)
                        }
                    }
                },
            );

            networks.push(data_net);
        }

        TensorNetworkSet { networks }
    }
}

#[cfg(feature = "shadowing")]
impl<
    T,
    S: TensorStructure + Clone,
    K: Clone,
    FK: Clone,
    Aind: AbsInd,
    Store: TensorScalarStore<Tensor = DataTensor<usize, S>, Scalar = usize> + Clone,
> EvalTensorNetworkSet<T, S, K, FK, Aind, Store>
{
    #[allow(clippy::type_complexity)]
    pub fn evaluate(
        &mut self,
        params: &[T],
    ) -> TensorNetworkSet<Store::Store<DataTensor<T, S>, T>, K, FK, Aind>
    where
        T: Real + SingleFloat,
        S: TensorStructure + Clone,
    {
        let zero = params[0].zero();
        let mut data = vec![zero; self.len];

        let mut networks = vec![];

        self.shared_data.evaluate(params, &mut data);

        for net in self.networks.iter() {
            let data_net = net.map_ref(
                |s| data[*s].clone(),
                |p| {
                    let structure = p.structure().clone();
                    match &p {
                        DataTensor::Dense(d) => {
                            let mut t_data = vec![];
                            for (_, &a) in d.flat_iter() {
                                t_data.push(data[a].clone());
                            }
                            DataTensor::Dense(DenseTensor::from_data(t_data, structure).unwrap())
                        }
                        DataTensor::Sparse(s) => {
                            let mut t = SparseTensor::empty(structure, T::new_zero());
                            for (i, &a) in s.flat_iter() {
                                t.set_flat(i, data[a].clone()).unwrap();
                            }
                            DataTensor::Sparse(t)
                        }
                    }
                },
            );

            networks.push(data_net);
        }

        TensorNetworkSet { networks }
    }

    /// Create a C++ code representation of the evaluation tree tensor.
    /// With `inline_asm` set to any value other than `None`,
    /// high-performance inline ASM code will be generated for most
    /// evaluation instructions. This often gives better performance than
    /// the `O3` optimization level and results in very fast compilation.
    #[allow(clippy::type_complexity)]
    pub fn export_cpp<F: CompiledNumber>(
        &self,
        path: impl AsRef<std::path::Path>,
        function_name: &str,
        settings: ExportSettings,
    ) -> Result<SharedTensorNetworkSet<ExportedCode<F>, S, K, FK, Aind, Store>, std::io::Error>
    where
        T: ExportNumber + SingleFloat,
    {
        Ok(SharedTensorNetworkSet {
            networks: self.networks.clone(),
            shared_data: self.shared_data.export_cpp(path, function_name, settings)?,
            len: self.len,
        })
    }
}

#[cfg(feature = "shadowing")]
impl<F: CompiledNumber, S: TensorStructure + Clone, K: Clone, FK: Clone, Aind: AbsInd>
    SharedTensorNetworkSet<ExportedCode<F>, S, K, FK, Aind>
{
    pub fn compile(
        &self,
        out: impl AsRef<std::path::Path>,
        options: CompileOptions,
    ) -> Result<SharedTensorNetworkSet<CompiledCode<F>, S, K, FK, Aind>, std::io::Error> {
        Ok(SharedTensorNetworkSet {
            networks: self.networks.clone(),
            shared_data: self.shared_data.compile(out, options)?,
            len: self.len,
        })
    }
}

#[cfg(feature = "shadowing")]
impl<F: CompiledNumber, S: TensorStructure + Clone, K: Clone, FK: Clone, Aind: AbsInd>
    SharedTensorNetworkSet<CompiledCode<F>, S, K, FK, Aind>
{
    pub fn load(&self) -> Result<SharedTensorNetworkSet<F::Evaluator, S, K, FK, Aind>, String> {
        Ok(SharedTensorNetworkSet {
            networks: self.networks.clone(),
            shared_data: self.shared_data.load()?,
            len: self.len,
        })
    }
}

#[cfg(feature = "shadowing")]
impl<
    S: TensorStructure + Clone,
    K: Clone,
    FK: Clone,
    Store: TensorScalarStore<Tensor = DataTensor<usize, S>, Scalar = usize> + Clone,
    Aind: AbsInd,
> CompiledTensorNetworkSet<S, K, FK, Aind, Store>
{
    #[allow(clippy::type_complexity, clippy::result_large_err)]
    pub fn evaluate(
        &mut self,
        params: &[crate::algebra::complex::Complex<f64>],
    ) -> TensorNetworkSet<
        Store::Store<
            DataTensor<crate::algebra::complex::Complex<f64>, S>,
            crate::algebra::complex::Complex<f64>,
        >,
        K,
        FK,
        Aind,
    >
    where
        S: TensorStructure + Clone,
    {
        let zero = crate::algebra::complex::Complex::default();
        let mut data = vec![zero; self.len];

        let mut networks = vec![];

        self.shared_data.evaluate(params, &mut data);
        for net in self.networks.iter() {
            let data_net = net.map_ref(
                |s| data[*s],
                |p| {
                    let structure = p.structure().clone();
                    match &p {
                        DataTensor::Dense(d) => {
                            let mut t_data = vec![];
                            for (_, &a) in d.flat_iter() {
                                t_data.push(data[a]);
                            }
                            DataTensor::Dense(DenseTensor::from_data(t_data, structure).unwrap())
                        }
                        DataTensor::Sparse(s) => {
                            let mut t = SparseTensor::empty(structure, Complex::new_zero());
                            for (i, &a) in s.flat_iter() {
                                t.set_flat(i, data[a]).unwrap();
                            }
                            DataTensor::Sparse(t)
                        }
                    }
                },
            );

            networks.push(data_net);
        }
        TensorNetworkSet { networks }
    }
}
