
use graph::{NAdd, NMul, NetworkEdge, NetworkGraph, NetworkLeaf, NetworkNode, NetworkOp};
use linnet::half_edge::NodeIndex;
use serde::{Deserialize, Serialize};

use library::{Library, LibraryError};

use crate::algebra::algebraic_traits::RefOne;
use crate::contraction::Contract;
use crate::network::library::{DummyKey, FunctionLibrary, FunctionLibraryError, LibraryTensor};
use crate::structure::abstract_index::AbstractIndex;
use crate::structure::permuted::PermuteTensor;
// use crate::shadowing::Concretize;
use crate::structure::representation::LibrarySlot;
use crate::structure::slot::{AbsInd, IsAbstractSlot};
use crate::structure::{HasName, PermutedStructure, StructureError};
use std::borrow::Cow;
use std::fmt::Display;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
use store::{NetworkStore, TensorScalarStore, TensorScalarStoreMapping};
use thiserror::Error;
// use log::trace;
use eyre::eyre;
use crate::{
    contraction::ContractionError,
    structure::{CastStructure, HasStructure, ScalarTensor, TensorStructure},
};

// use eyre::Result;

use std::{convert::Infallible, fmt::Debug};

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    PartialEq,
    Eq,
)]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
pub struct Network<S, LibKey, FunKey, Aind = AbstractIndex> {
    pub graph: NetworkGraph<LibKey, FunKey, Aind>,
    pub store: S,
    pub state: NetworkState,
}

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    PartialEq,
    Eq,
    Copy,
)]
pub enum NetworkState {
    PureScalar,
    Tensor,
    SelfDualTensor,
    Scalar,
}

impl NetworkState {
    pub fn is_scalar(&self) -> bool {
        matches!(self, NetworkState::Scalar | NetworkState::PureScalar)
    }

    pub fn is_tensor(&self) -> bool {
        matches!(self, NetworkState::Tensor | NetworkState::SelfDualTensor)
    }

    pub fn pow(self, pow: i8) -> Self {
        match self {
            NetworkState::PureScalar => NetworkState::PureScalar,
            NetworkState::Scalar => NetworkState::Scalar,
            NetworkState::SelfDualTensor => {
                if pow % 2 == 0 {
                    NetworkState::Scalar
                } else {
                    NetworkState::SelfDualTensor
                }
            }
            NetworkState::Tensor => panic!("Cannot have integer power of non-self dual tensor"),
        }
    }

    pub fn is_compatible(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (NetworkState::Tensor, NetworkState::Tensor)
                | (NetworkState::SelfDualTensor, NetworkState::SelfDualTensor)
                | (NetworkState::Scalar, NetworkState::Scalar)
                | (NetworkState::PureScalar, NetworkState::Scalar)
                | (NetworkState::Scalar, NetworkState::PureScalar)
                | (NetworkState::PureScalar, NetworkState::PureScalar)
        )
    }
}
impl MulAssign for NetworkState {
    fn mul_assign(&mut self, rhs: Self) {
        // println!("{self:?} *={rhs:?}");
        *self = match (*self, rhs) {
            (NetworkState::PureScalar, NetworkState::PureScalar) => NetworkState::PureScalar,
            (NetworkState::PureScalar, NetworkState::Scalar) => NetworkState::Scalar,
            (NetworkState::Tensor, _) => NetworkState::Tensor,
            (NetworkState::Scalar, NetworkState::PureScalar) => NetworkState::Scalar,
            (_, NetworkState::Tensor) => NetworkState::Tensor,
            (NetworkState::SelfDualTensor, _) => NetworkState::SelfDualTensor,
            (_, NetworkState::SelfDualTensor) => NetworkState::SelfDualTensor,
            (NetworkState::Scalar, NetworkState::Scalar) => NetworkState::Scalar,
        }
    }
}

impl AddAssign for NetworkState {
    fn add_assign(&mut self, rhs: Self) {
        // println!("{self:?} *={rhs:?}");
        *self = match (*self, rhs) {
            (NetworkState::PureScalar, NetworkState::PureScalar) => NetworkState::PureScalar,
            (NetworkState::PureScalar, NetworkState::Scalar) => NetworkState::Scalar,
            (a, b) => {
                assert_eq!(a, b, "Cannot add incompatible network states:{a:?} + {b:?}");
                a
            }
        }
    }
}

// pub type TensorNetwork<T, S, Str: TensorScalarStore<Tensor = T, Scalar = S>, K> = Network<Str, K>;

// pub struct TensorNetwork<
//     T,
//     S,
//     K,
//     Str: TensorScalarStore<Tensor = T, Scalar = S> = NetworkStore<T, S>,
// > {
//     net: Network<Str, K>,
// }

pub mod graph;
pub mod library;
pub mod set;
pub mod store;

impl<S: TensorScalarStoreMapping, K: Clone, FK: Clone, Aind: AbsInd> TensorScalarStoreMapping
    for Network<S, K, FK, Aind>
{
    type Store<U, V> = Network<S::Store<U, V>, K, FK, Aind>;
    type Scalar = S::Scalar;
    type Tensor = S::Tensor;

    fn iter_scalars(&self) -> impl Iterator<Item = &Self::Scalar> {
        self.store.iter_scalars()
    }

    fn iter_tensors(&self) -> impl Iterator<Item = &Self::Tensor> {
        self.store.iter_tensors()
    }

    fn iter_scalars_mut(&mut self) -> impl Iterator<Item = &mut Self::Scalar> {
        self.store.iter_scalars_mut()
    }
    fn iter_tensors_mut(&mut self) -> impl Iterator<Item = &mut Self::Tensor> {
        self.store.iter_tensors_mut()
    }

    fn map<U, V>(
        self,
        scalar_map: impl FnMut(Self::Scalar) -> U,
        tensor_map: impl FnMut(Self::Tensor) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map(scalar_map, tensor_map),
            graph: self.graph,
            state: self.state,
        }
    }

    fn map_result<U, V, Er>(
        self,
        scalar_map: impl FnMut(Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self.store.map_result(scalar_map, tensor_map)?,
            graph: self.graph,
            state: self.state,
        })
    }

    fn map_ref<'a, U, V>(
        &'a self,
        scalar_map: impl FnMut(&'a Self::Scalar) -> U,
        tensor_map: impl FnMut(&'a Self::Tensor) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map_ref(scalar_map, tensor_map),
            graph: self.graph.clone(),
            state: self.state,
        }
    }

    fn map_ref_result<U, V, Er>(
        &self,
        scalar_map: impl FnMut(&Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(&Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self.store.map_ref_result(scalar_map, tensor_map)?,
            graph: self.graph.clone(),
            state: self.state,
        })
    }

    fn map_ref_enumerate<U, V>(
        &self,
        scalar_map: impl FnMut((usize, &Self::Scalar)) -> U,
        tensor_map: impl FnMut((usize, &Self::Tensor)) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map_ref_enumerate(scalar_map, tensor_map),
            graph: self.graph.clone(),
            state: self.state,
        }
    }

    fn map_ref_result_enumerate<U, V, Er>(
        &self,
        scalar_map: impl FnMut((usize, &Self::Scalar)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &Self::Tensor)) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self
                .store
                .map_ref_result_enumerate(scalar_map, tensor_map)?,
            graph: self.graph.clone(),
            state: self.state,
        })
    }

    fn map_ref_mut<U, V>(
        &mut self,
        scalar_map: impl FnMut(&mut Self::Scalar) -> U,
        tensor_map: impl FnMut(&mut Self::Tensor) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map_ref_mut(scalar_map, tensor_map),
            graph: self.graph.clone(),
            state: self.state,
        }
    }

    fn map_ref_mut_result<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut(&mut Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(&mut Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self.store.map_ref_mut_result(scalar_map, tensor_map)?,
            graph: self.graph.clone(),
            state: self.state,
        })
    }

    fn map_ref_mut_enumerate<U, V>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut Self::Scalar)) -> U,
        tensor_map: impl FnMut((usize, &mut Self::Tensor)) -> V,
    ) -> Self::Store<V, U> {
        Network {
            store: self.store.map_ref_mut_enumerate(scalar_map, tensor_map),
            graph: self.graph.clone(),
            state: self.state,
        }
    }

    fn map_ref_mut_result_enumerate<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut Self::Scalar)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &mut Self::Tensor)) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er> {
        Ok(Network {
            store: self
                .store
                .map_ref_mut_result_enumerate(scalar_map, tensor_map)?,
            graph: self.graph.clone(),
            state: self.state,
        })
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> Default for Network<S, K, FK, Aind> {
    fn default() -> Self {
        Self::one()
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> NMul for Network<S, K, FK, Aind> {
    type Output = Self;
    fn n_mul<I: IntoIterator<Item = Self>>(self, iter: I) -> Self::Output {
        let mut store = self.store;
        let mut state = self.state;

        let items = iter.into_iter().map(|mut a| {
            a.graph.shift_scalars(store.n_scalars());
            a.graph.shift_tensors(store.n_tensors());
            store.extend(a.store);

            state *= a.state;
            a.graph
        });

        let graph = self.graph.n_mul(items);

        if state.is_tensor() && graph.dangling_indices().is_empty() {
            state = NetworkState::Scalar;
        }

        Network {
            graph,
            store,
            state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> Mul for Network<S, K, FK, Aind> {
    type Output = Self;
    fn mul(mut self, mut other: Self) -> Self::Output {
        let mut store = self.store;

        other.graph.shift_scalars(store.n_scalars());
        other.graph.shift_tensors(store.n_tensors());
        store.extend(other.store);
        self.state *= other.state;

        Network {
            graph: self.graph * other.graph,
            store,
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> MulAssign
    for Network<S, K, FK, Aind>
{
    fn mul_assign(&mut self, mut rhs: Self) {
        rhs.graph.shift_scalars(self.store.n_scalars());
        rhs.graph.shift_tensors(self.store.n_tensors());
        self.store.extend(rhs.store);
        self.state *= rhs.state;
        self.graph *= rhs.graph;
    }
}

impl<T: TensorStructure, S, FK: Debug, K: Debug, Aind: AbsInd> MulAssign<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    fn mul_assign(&mut self, rhs: T) {
        *self *= Network::from_tensor(rhs);
    }
}

impl<T: TensorStructure, S, FK: Debug, K: Debug, Aind: AbsInd> Mul<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    type Output = Self;
    fn mul(self, other: T) -> Self::Output {
        let mut store = self.store;

        let mut other = Network::from_tensor(other);

        other.graph.shift_scalars(store.n_scalars());
        other.graph.shift_tensors(store.n_tensors());
        store.extend(other.store);

        Network {
            graph: self.graph * other.graph,
            store,
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> Add for Network<S, K, FK, Aind> {
    type Output = Self;
    fn add(self, mut other: Self) -> Self::Output {
        let mut store = self.store;

        other.graph.shift_scalars(store.n_scalars());
        other.graph.shift_tensors(store.n_tensors());
        store.extend(other.store);

        Network {
            graph: self.graph + other.graph,
            store,
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> AddAssign
    for Network<S, K, FK, Aind>
{
    fn add_assign(&mut self, mut rhs: Self) {
        rhs.graph.shift_scalars(self.store.n_scalars());
        rhs.graph.shift_tensors(self.store.n_tensors());
        self.store.extend(rhs.store);

        self.graph += rhs.graph;
    }
}

impl<T: TensorStructure, S, FK: Debug, K: Debug, Aind: AbsInd> AddAssign<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    fn add_assign(&mut self, rhs: T) {
        *self += Network::from_tensor(rhs);
    }
}

impl<T: TensorStructure, S, FK: Debug, K: Debug, Aind: AbsInd> Add<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    type Output = Self;
    fn add(mut self, other: T) -> Self::Output {
        self += other;
        self
    }
}

impl<T: TensorStructure, FK: Debug, K: Debug, Aind: AbsInd> Add<i8>
    for Network<NetworkStore<T, i8>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    type Output = Self;
    fn add(mut self, other: i8) -> Self::Output {
        let mut other = Network::from_scalar(other);
        other.graph.shift_tensors(self.store.n_tensors());
        other.graph.shift_tensors(self.store.n_tensors());

        self.store.extend(other.store);
        Network {
            graph: self.graph + other.graph,
            store: self.store,
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> NAdd for Network<S, K, FK, Aind> {
    type Output = Self;
    fn n_add<I: IntoIterator<Item = Self>>(self, iter: I) -> Self::Output {
        let mut store = self.store;

        let mut state = self.state;

        let items = iter.into_iter().map(|mut a| {
            a.graph.shift_scalars(store.n_scalars());
            a.graph.shift_tensors(store.n_tensors());
            store.extend(a.store);

            state += a.state;
            a.graph
        });

        Network {
            graph: self.graph.n_add(items),
            store,
            state,
        }
    }
}

impl<S: TensorScalarStore, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Neg
    for Network<S, K, FK, Aind>
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            store: self.store,
            graph: self.graph.neg(),
            state: self.state,
        }
    }
}

impl<S: TensorScalarStore, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Sub
    for Network<S, K, FK, Aind>
{
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<S: TensorScalarStore, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> SubAssign
    for Network<S, K, FK, Aind>
{
    fn sub_assign(&mut self, mut rhs: Self) {
        rhs.graph.shift_scalars(self.store.n_scalars());
        rhs.graph.shift_tensors(self.store.n_tensors());
        self.store.extend(rhs.store);

        self.graph -= rhs.graph
    }
}

impl<T: TensorStructure, S, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> SubAssign<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    fn sub_assign(&mut self, rhs: T) {
        *self -= Network::from_tensor(rhs)
    }
}

impl<T: TensorStructure, S, K: Clone + Debug, FK: Clone + Debug, Aind: AbsInd> Sub<T>
    for Network<NetworkStore<T, S>, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    type Output = Self;
    fn sub(mut self, other: T) -> Self::Output {
        self -= other;
        self
    }
}

impl<S: TensorScalarStore, FK: Debug, K: Debug, Aind: AbsInd> Network<S, K, FK, Aind> {
    pub fn pow(self, pow: i8) -> Self {
        Self {
            store: self.store,
            graph: self.graph.pow(pow),
            state: self.state.pow(pow),
        }
    }
    pub fn fun(self, key: FK) -> Self {
        let graph = self.graph.function(key);
        Self {
            store: self.store,
            state: graph.state(),
            graph,
        }
    }

    pub fn from_scalar(scalar: S::Scalar) -> Self {
        let mut store = S::default();
        let id = store.add_scalar(scalar);
        Network {
            graph: NetworkGraph::scalar(id),
            store,

            state: NetworkState::PureScalar,
        }
    }

    pub fn merge_ops(&mut self)
    where
        K: Clone,
    {
        self.graph.merge_ops();
    }

    pub fn from_tensor(tensor: S::Tensor) -> Self
    where
        S::Tensor: TensorStructure,
        <S::Tensor as TensorStructure>::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let mut store = S::default();

        let state = if tensor.is_scalar() {
            NetworkState::Scalar
        } else if tensor.is_fully_self_dual() {
            // tensor.structure().dual();
            NetworkState::SelfDualTensor
        } else {
            NetworkState::Tensor
        };
        let id = store.add_tensor(tensor);

        Network {
            graph: NetworkGraph::tensor(store.get_tensor(id), NetworkLeaf::LocalTensor(id)),
            store,
            state,
        }
    }

    pub fn library_tensor<T>(tensor: &T, key: PermutedStructure<K>) -> Self
    where
        T: TensorStructure,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let state = if tensor.is_scalar() {
            NetworkState::Scalar
        } else if tensor.is_fully_self_dual() {
            NetworkState::SelfDualTensor
        } else {
            NetworkState::Tensor
        };
        Network {
            graph: NetworkGraph::tensor(tensor, NetworkLeaf::LibraryKey(key)),
            store: S::default(),
            state,
        }
    }

    pub fn one() -> Self {
        Network {
            graph: NetworkGraph::one(),
            store: S::default(),
            state: NetworkState::PureScalar,
        }
    }

    pub fn zero() -> Self {
        Network {
            graph: NetworkGraph::zero(),
            store: S::default(),
            state: NetworkState::PureScalar,
        }
    }
}

#[derive(Error, Debug)]
pub enum TensorNetworkError<K: Display, FK: Display> {
    #[error("Slot edge to prod node")]
    SlotEdgeToProdNode,
    #[error("Slot edge to scalar node")]
    SlotEdgeToScalarNode,
    #[error("More than one neg")]
    MoreThanOneNeg,
    #[error("Childless neg")]
    ChildlessNeg,
    #[error("Contraction Error:{0}")]
    ContractionError(#[from] ContractionError),
    #[error("Scalar connected by a slot edge")]
    ScalarSlotEdge,
    #[error("Structure Error:{0}")]
    StructErr(#[from] StructureError),
    #[error("LibraryError:{0}")]
    LibErr(#[from] LibraryError<K>),
    #[error("FunctionLibraryError:{0}")]
    FunLibErr(#[from] FunctionLibraryError<FK>),
    #[error("Non tensor node still present")]
    NonTensorNodePresent,
    #[error("Negative non-even power on non-scalar node:{0}")]
    NegativeExponentNonScalar(String),
    #[error("Too many arguments for function:{0}")]
    TooManyArgsFunction(String),
    #[error("Non self-dual tensor power{0}")]
    InvalidDotFunction(String),
    #[error("Invalid dot function{0}")]
    NonSelfDualTensorPower(String),
    #[error("invalid resulting node{0}")]
    InvalidResultNode(NetworkNode<DummyKey, FK>),
    #[error("internal edge still present, contract it first")]
    InternalEdgePresent,
    #[error("uncontracted scalar")]
    UncontractedScalar,
    #[error("Cannot contract edge between {0} and {1}")]
    CannotContractEdgeBetween(NetworkNode<K, FK>, NetworkNode<K, FK>),
    #[error("no nodes in the graph")]
    NoNodes,
    #[error("no scalar present")]
    NoScalar,
    #[error("more than one node in the graph")]
    MoreThanOneNode,
    #[error("is not scalar output")]
    NotScalarOutput,
    #[error("failed scalar multiplication")]
    FailedScalarMul,
    #[error("scalar field is empty")]
    ScalarFieldEmpty,
    #[error("not all scalars: {0}")]
    NotAllScalars(String),
    #[error("try to sum scalar with library tensor: {0}")]
    ScalarLibSum(String),
    #[error("try to sum scalar with a tensor: {0}")]
    SumScalarTensor(String),
    #[error("Incompatible summands: {0}")]
    IncompatibleSummand(String),
    #[error("failed to contract")]
    FailedContract(ContractionError),
    #[error("negative exponent not yet supported")]
    NegativeExponent,
    #[error("failed to contract: {0}")]
    FailedContractMsg(String),
    #[error(transparent)]
    Other(#[from] eyre::Error),
    #[error("Io error")]
    InOut(#[from] std::io::Error),
    #[error("Infallible")]
    Infallible,
}

impl<K: Display, FK: Display> From<Infallible> for TensorNetworkError<K, FK> {
    fn from(_: Infallible) -> Self {
        TensorNetworkError::Infallible
    }
}

pub enum TensorOrScalarOrKey<T, S, K, Aind> {
    Tensor {
        tensor: T,
        graph_slots: Vec<LibrarySlot<Aind>>,
    },
    Scalar(S),
    Key {
        key: K,
        nodeid: NodeIndex,
    },
}

pub enum ExecutionResult<T> {
    One,
    Zero,
    Val(T),
}

impl<T: Display> Display for ExecutionResult<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExecutionResult::One => write!(f, "One"),
            ExecutionResult::Zero => write!(f, "Zero"),
            ExecutionResult::Val(val) => write!(f, "{}", val),
        }
    }
}

impl<
    T: TensorStructure,
    S,
    K: Display + Debug,
    FK: Display + Debug,
    Str: TensorScalarStore<Tensor = T, Scalar = S>,
    Aind: AbsInd,
> Network<Str, K, FK, Aind>
where
    T::Slot: IsAbstractSlot<Aind = Aind>,
{
    pub fn validate(&self)
    where
        K: TensorStructure,
    {
        for (n, _neigh, v) in self.graph.graph.iter_nodes() {
            match v {
                NetworkNode::Leaf(NetworkLeaf::LibraryKey(k)) => {
                    let reps = self
                        .graph
                        .slots(n)
                        .into_iter()
                        .map(|s| s.rep())
                        .collect::<Vec<_>>();
                    // let p = Permutation::sort(&reps);

                    let n_reps = k
                        .structure
                        .external_reps_iter()
                        .map(|r| r.to_lib())
                        .collect::<Vec<_>>();
                    // let q = Permutation::sort(&n_reps);
                    // println!("p{p}q{q}");
                    assert_eq!(n_reps, reps);
                }
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(k)) => {
                    let reps = self
                        .graph
                        .slots(n)
                        .into_iter()
                        .map(|s| s.rep())
                        .collect::<Vec<_>>();
                    let n_reps = self
                        .store
                        .get_tensor(*k)
                        .external_reps_iter()
                        .map(|r| r.to_lib())
                        .collect::<Vec<_>>();
                    assert_eq!(n_reps, reps);
                }
                _ => {}
            }
        }
    }

    #[allow(clippy::result_large_err, clippy::type_complexity)]
    pub fn result(
        &self,
    ) -> Result<
        ExecutionResult<TensorOrScalarOrKey<&T, &S, &PermutedStructure<K>, Aind>>,
        TensorNetworkError<K, FK>,
    >
    where
        FK: Clone,
    {
        let (node, nid, graph_slots) = self.graph.result()?;

        match node {
            NetworkNode::Leaf(l) => match l {
                NetworkLeaf::LibraryKey(k) => Ok(ExecutionResult::Val(TensorOrScalarOrKey::Key {
                    key: k,
                    nodeid: nid,
                })),
                NetworkLeaf::LocalTensor(t) => {
                    Ok(ExecutionResult::Val(TensorOrScalarOrKey::Tensor {
                        tensor: self.store.get_tensor(*t),
                        graph_slots,
                    }))
                }
                NetworkLeaf::Scalar(t) => Ok(ExecutionResult::Val(TensorOrScalarOrKey::Scalar(
                    self.store.get_scalar(*t),
                ))),
            },
            NetworkNode::Op(o) => match o {
                NetworkOp::Product => Ok(ExecutionResult::One),
                NetworkOp::Sum => Ok(ExecutionResult::Zero),
                o => Err(TensorNetworkError::InvalidResultNode(NetworkNode::Op(
                    o.clone(),
                ))),
            },
        }
    }

    #[allow(clippy::result_large_err)]
    pub fn result_tensor<'a, LT, L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>>(
        &'a self,
        lib: &L,
    ) -> Result<ExecutionResult<Cow<'a, T>>, TensorNetworkError<K, FK>>
    where
        S: 'a,
        T: Clone + ScalarTensor + HasStructure,
        K: Display + Debug,
        FK: Display + Debug + Clone,
        LT: TensorStructure<Indexed = T> + Clone + LibraryTensor<WithIndices = T>,
        T: PermuteTensor<Permuted = T>,
        for<'b> &'b S: Into<T::Scalar>,
        <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
            IsAbstractSlot<Aind = Aind>,
    {
        Ok(match self.result()? {
            ExecutionResult::One => ExecutionResult::One,
            ExecutionResult::Zero => ExecutionResult::Zero,
            ExecutionResult::Val(v) => ExecutionResult::Val(match v {
                TensorOrScalarOrKey::Tensor { tensor, .. } => Cow::Borrowed(tensor),
                TensorOrScalarOrKey::Scalar(s) => Cow::Owned(T::new_scalar(s.into())),
                TensorOrScalarOrKey::Key { nodeid, .. } => {
                    let less = self.graph.get_lib_data(lib, nodeid).unwrap();

                    Cow::Owned(less)
                }
            }),
        })
    }

    #[allow(clippy::result_large_err)]
    pub fn result_scalar<'a>(
        &'a self,
    ) -> Result<ExecutionResult<Cow<'a, S>>, TensorNetworkError<K, FK>>
    where
        T: Clone + ScalarTensor + 'a,
        T::Scalar: Into<S>,
        K: Display,
        FK: Display + Clone,
        S: Clone,
    {
        Ok(match self.result()? {
            ExecutionResult::One => ExecutionResult::One,
            ExecutionResult::Zero => ExecutionResult::Zero,
            ExecutionResult::Val(v) => ExecutionResult::Val(match v {
                TensorOrScalarOrKey::Tensor { tensor: t, .. } => Cow::Owned(
                    t.clone()
                        .scalar()
                        .ok_or(TensorNetworkError::NoScalar)?
                        .into(),
                ),
                TensorOrScalarOrKey::Scalar(s) => Cow::Borrowed(s),
                TensorOrScalarOrKey::Key { .. } => return Err(TensorNetworkError::NoScalar),
            }),
        })
    }

    pub fn cast<U>(self) -> Network<Str::Store<U, S>, K, FK, Aind>
    where
        K: Clone,
        FK: Clone,
        T: CastStructure<U> + HasStructure,
        T::Structure: TensorStructure,
        U: HasStructure,
        U::Structure: From<T::Structure> + TensorStructure<Slot = T::Slot>,
    {
        self.map(|a| a, |t| t.cast_structure())
    }
}

pub trait StructureLessDisplay {
    fn display(&self) -> String {
        String::new()
    }
}

impl<S: HasName> StructureLessDisplay for S
where
    S::Args: StructureLessDisplay,
    S::Name: Display,
{
    fn display(&self) -> String {
        format!(
            "{}({})",
            self.name().map(|t| t.to_string()).unwrap_or_default(),
            self.args().map(|t| t.display()).unwrap_or_default()
        )
    }
}

impl<S, K: Display + Debug, FK: Display + Debug, Aind: AbsInd> Network<S, K, FK, Aind> {
    pub fn dot(&self) -> std::string::String {
        self.graph.dot()
    }

    pub fn dot_pretty(&self) -> std::string::String
    where
        S: TensorScalarStore,
        S::Scalar: Display,
        K: StructureLessDisplay,
        S::Tensor: StructureLessDisplay,
    {
        self.graph.dot_impl(
            |i| {
                let ss = &self.store.get_scalar(i);
                format!("{}:{}", i, ss)
            },
            |k| k.display(),
            |t| {
                let tt = &self.store.get_tensor(t);
                tt.display()
            },
            |fk| fk.to_string(),
        )
    }
}

impl<T, S, FK: Debug, K: Debug, Aind: AbsInd> Network<NetworkStore<T, S>, K, FK, Aind> {
    pub fn dot_display_impl(
        &self,
        scalar_disp: impl Fn(&S) -> String,
        library_disp: impl Fn(&K) -> Option<String>,
        tensor_disp: impl Fn(&T) -> String,
        function_disp: impl Fn(&FK) -> String,
    ) -> std::string::String {
        self.graph.graph.dot_impl(
            &self.graph.graph.full_filter(),
            "",
            &|_| None,
            &|e| {
                if let NetworkEdge::Slot(s) = e {
                    Some(format!("label=\"{s}\""))
                } else {
                    None
                }
            },
            &|n| match n {
                NetworkNode::Leaf(l) => match l {
                    NetworkLeaf::LibraryKey(l) => {
                        // if let Ok(v) = lib.get(l) {
                        Some(format!("label = \"L:{}\"", library_disp(&l.structure)?))
                        // } else {
                        // None
                        // }
                    }
                    NetworkLeaf::LocalTensor(l) => Some(format!(
                        "label = \"T:{}\"",
                        tensor_disp(self.store.get_tensor(*l))
                    )),
                    NetworkLeaf::Scalar(s) => Some(format!(
                        "label = \"S:{}\"",
                        scalar_disp(self.store.get_scalar(*s))
                    )),
                },
                NetworkNode::Op(o) => {
                    Some(format!("label = \"{}\"", o.display_with(&function_disp)))
                }
            },
        )
        // self.graph.dot()
    }
}

// use log::trace;
#[cfg(feature = "shadowing")]
pub mod parsing;
// use log::trace;
pub mod contract;
pub use contract::{
    ContractScalars, ContractionStrategy, SingleSmallestDegree, SmallestDegree, SmallestDegreeIter,
};
pub trait ExecutionStrategy<E, FL, L, K, FK, Aind>
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
{
    /// Run the entire contraction to one leaf.
    #[allow(clippy::result_large_err)]
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;
}

pub struct Sequential;

pub struct Steps<const N: usize> {}
pub struct StepsDebug<const N: usize> {}

impl<const N: usize, E, L, FL, FK: Debug, K: Debug, Aind: AbsInd>
    ExecutionStrategy<E, FL, L, K, FK, Aind> for StepsDebug<N>
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
    K: Clone,
    FK: Clone,
{
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        for _ in 0..N {
            // find the *one* ready op
            if let Some((extracted_graph, op)) = graph.extract_next_ready_op() {
                println!(
                    "Extracted_graph: {}",
                    extracted_graph.dot_impl(
                        |s| s.to_string(),
                        |_| "".to_string(),
                        |s| s.to_string(),
                        |_| "".to_string()
                    )
                );
                println!(
                    "Graph: {}",
                    graph.dot_impl(
                        |s| s.to_string(),
                        |_| "".to_string(),
                        |s| s.to_string(),
                        |_| "".to_string()
                    )
                );
                // execute + splice
                let replacement = executor.execute::<C>(extracted_graph, lib, fnlib, op)?;
                println!(
                    "Replacement Graph: {}",
                    replacement.dot_impl(
                        |s| s.to_string(),
                        |_| "".to_string(),
                        |s| s.to_string(),
                        |_| "".to_string()
                    )
                );

                graph.splice_descendents_of(replacement);
            }
        }

        Ok(())
    }
}

impl<const N: usize, E, FL, L, K, FK, Aind: AbsInd> ExecutionStrategy<E, FL, L, K, FK, Aind>
    for Steps<N>
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
    K: Clone + Debug,
    FK: Clone + Debug,
{
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        for _ in 0..N {
            // find the *one* ready op
            if let Some((extracted_graph, op)) = graph.extract_next_ready_op() {
                // execute + splice
                let replacement = executor.execute::<C>(extracted_graph, lib, fnlib, op)?;
                graph.splice_descendents_of(replacement);
            }
        }

        Ok(())
    }
}

impl<E, L, FL, K, FK, Aind: AbsInd> ExecutionStrategy<E, FL, L, K, FK, Aind> for Sequential
where
    E: ExecuteOp<FL, L, K, FK, Aind>,
    FK: Clone + Debug,
    K: Clone + Debug,
{
    fn execute_all<C: ContractionStrategy<E, L, K, FK, Aind>>(
        executor: &mut E,
        graph: &mut NetworkGraph<K, FK, Aind>,
        lib: &L,
        fnlib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display,
    {
        while {
            // find the *one* ready op
            if let Some((extracted_graph, op)) = graph.extract_next_ready_op() {
                // execute + splice
                let replacement = executor.execute::<C>(extracted_graph, lib, fnlib, op)?;
                graph.splice_descendents_of(replacement);
                true
            } else {
                false
            }
        } {}

        Ok(())
    }
}

// 2b) Parallel: batch‐execute all ready ops, then splice serially.
// pub struct Parallel;
// impl<E, K> ExecutionStrategy<E, K> for Parallel
// where
//     E: ExecuteOp<K> + Clone + Send + Sync,
//     K: Clone + Send + Sync,
// {
//     fn contract_all<L: Library<Key = K> + Sync>(
//         &self,
//         executor: &mut E,
//         graph: &mut NetworkGraph<K>,
//         lib: &L,
//     ) {
//         loop {
//             // 1) collect *all* ready ops this round
//             let ready = graph.find_all_ready_ops();

//             if ready.is_empty() {
//                 break;
//             }

//             // 2) execute them in parallel
//             let results: Vec<(NodeIndex, NetworkGraph<K>)> = ready
//                 .into_par_iter()
//                 .map(|(nid, op, leaves)| {
//                     let mut local = executor.clone();
//                     let replacement = local.execute(lib, op, &leaves);
//                     (nid, replacement)
//                 })
//                 .collect();

//             // 3) splice back sequentially
//             for (nid, replacement) in results {
//                 graph.splice_descendents_of(nid, replacement);
//             }
//         }
//     }
// }

pub trait ExecuteOp<FL, L, K, FK, Aind>: Sized {
    // type LibStruct;
    #[allow(clippy::result_large_err)]
    fn execute<C: ContractionStrategy<Self, L, K, FK, Aind>>(
        &mut self,
        graph: NetworkGraph<K, FK, Aind>,
        lib: &L,
        fn_lib: &FL,
        op: NetworkOp<FK>,
    ) -> Result<NetworkGraph<K, FK, Aind>, TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;
}

impl<S, Store: TensorScalarStore, K, FK, Aind: AbsInd> Network<Store, K, FK, Aind>
where
    Store::Tensor: HasStructure<Structure = S>,
{
    #[allow(clippy::result_large_err)]
    pub fn execute<
        Strat: ExecutionStrategy<Store, FL, L, K, FK, Aind>,
        C: ContractionStrategy<Store, L, K, FK, Aind>,
        LT,
        L,
        FL,
    >(
        &mut self,
        lib: &L,
        fn_lib: &FL,
    ) -> Result<(), TensorNetworkError<K, FK>>
    where
        K: Display + Clone + Debug,
        FK: Display + Clone + Debug,
        L: Library<S, Key = K, Value = PermutedStructure<LT>> + Sync,
        FL: FunctionLibrary<Store::Tensor, Store::Scalar, Key = FK>,
        LT: LibraryTensor<WithIndices = Store::Tensor>,
        Store: ExecuteOp<FL, L, K, FK, Aind>,
    {
        self.merge_ops();
        // println!("Hi");
        // println!("{}", self.graph.dot());
        // Ok(())
        Strat::execute_all::<C>(&mut self.store, &mut self.graph, lib, fn_lib)
    }
}

impl<
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Neg<Output = T>
        + Clone
        + Ref
        + Contract<LCM = T>
        + for<'a> AddAssign<T::Ref<'a>>
        + for<'a> AddAssign<LT::WithIndices>
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: Neg<Output = Sc>
        + RefOne
        + Div<Output = Sc>
        + for<'a> AddAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> AddAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref
        + for<'a> MulAssign<Sc::Ref<'a>>,
    K: Display + Debug + Clone,
    FK: Display + Debug,
    FL: FunctionLibrary<T, Sc, Key = FK>,
    Aind: AbsInd,
> ExecuteOp<FL, L, K, FK, Aind> for NetworkStore<T, Sc>
where
    LT::WithIndices: PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn execute<C: ContractionStrategy<Self, L, K, FK, Aind>>(
        &mut self,
        mut graph: NetworkGraph<K, FK, Aind>,
        lib: &L,
        fn_lib: &FL,
        op: NetworkOp<FK>,
    ) -> Result<NetworkGraph<K, FK, Aind>, TensorNetworkError<K, FK>> {
        graph.sync_order();
        match op {
            NetworkOp::Neg => {
                let ops = graph
                    .graph
                    .iter_nodes()
                    .find(|(_, _, d)| matches!(d, NetworkNode::Op(NetworkOp::Neg)));

                let (opid, children, _) = ops.unwrap();

                let mut child = None;
                for c in children {
                    if let Some(id) = graph.graph.involved_node_id(c)
                        && let NetworkNode::Leaf(l) = &graph.graph[id]
                    {
                        if child.is_some() {
                            return Err(TensorNetworkError::MoreThanOneNeg);
                        } else {
                            child = Some((id, l));
                        }
                    }
                }
                if let Some((child_id, leaf)) = child {
                    let new_node = match leaf {
                        NetworkLeaf::Scalar(s) => {
                            let s = self.scalar[*s].clone().neg();
                            let pos = self.scalar.len();
                            self.scalar.push(s);

                            NetworkLeaf::Scalar(pos)
                        }
                        NetworkLeaf::LibraryKey(_) => {
                            let inds = graph.get_lib_data(lib, child_id).unwrap();

                            let t = T::from(inds).neg();
                            let pos = self.tensors.len();
                            self.tensors.push(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        NetworkLeaf::LocalTensor(t) => {
                            let t = self.tensors[*t].clone().neg();
                            let pos = self.tensors.len();
                            self.tensors.push(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                    };
                    graph.identify_nodes_without_self_edges(
                        &[child_id, opid],
                        NetworkNode::Leaf(new_node),
                    );
                    Ok(graph)
                } else {
                    Err(TensorNetworkError::ChildlessNeg)
                }
            }
            NetworkOp::Product => {
                // println!("Doing Product");
                let (graph, _) = C::contract(self, graph, lib)?;
                Ok(graph)
            }
            NetworkOp::Sum => {
                // let mut op = None;
                let mut targets = Vec::new();
                let mut all_nodes = Vec::new();
                for (n, _, v) in graph.graph.iter_nodes() {
                    all_nodes.push(n);
                    if let NetworkNode::Leaf(l) = &v {
                        targets.push((n, l));
                    }
                }

                let (nf, first) = &targets[0];

                let new_node = match first {
                    NetworkLeaf::Scalar(s) => {
                        let mut accumulator = self.scalar[*s].clone();

                        for (_, t) in &targets[1..] {
                            match t {
                                NetworkLeaf::Scalar(s) => {
                                    accumulator += self.scalar[*s].refer();
                                }
                                NetworkLeaf::LocalTensor(t) => {
                                    if let Some(s) = self.tensors[*t].scalar_ref() {
                                        accumulator += s;
                                    } else {
                                        return Err(TensorNetworkError::NotAllScalars(
                                            "".to_string(),
                                        ));
                                    }
                                }
                                NetworkLeaf::LibraryKey { .. } => {
                                    return Err(TensorNetworkError::ScalarLibSum("".to_string()));
                                }
                            }
                        }

                        let pos = self.scalar.len();
                        self.scalar.push(accumulator);
                        NetworkLeaf::Scalar(pos)
                    }
                    NetworkLeaf::LocalTensor(t) => {
                        let mut accumulator = self.tensors[*t].clone();
                        if accumulator.is_scalar() {
                            let mut accumulator = Sc::from(accumulator.scalar().unwrap());

                            for (_, t) in &targets[1..] {
                                match t {
                                    NetworkLeaf::Scalar(s) => {
                                        accumulator += self.scalar[*s].refer();
                                    }
                                    NetworkLeaf::LocalTensor(t) => {
                                        if let Some(s) = self.tensors[*t].scalar_ref() {
                                            accumulator += s;
                                        } else {
                                            return Err(TensorNetworkError::NotAllScalars(
                                                "".to_string(),
                                            ));
                                        }
                                    }
                                    NetworkLeaf::LibraryKey { .. } => {
                                        return Err(TensorNetworkError::ScalarLibSum(
                                            "".to_string(),
                                        ));
                                    }
                                }
                            }

                            let pos = self.scalar.len();
                            self.scalar.push(accumulator);
                            NetworkLeaf::Scalar(pos)
                        } else {
                            for (nid, t) in &targets[1..] {
                                match t {
                                    NetworkLeaf::Scalar(_) => {
                                        return Err(TensorNetworkError::SumScalarTensor(
                                            "".to_string(),
                                        ));
                                    }
                                    NetworkLeaf::LocalTensor(t) => {
                                        accumulator += self.tensors[*t].refer();
                                    }
                                    NetworkLeaf::LibraryKey(_) => {
                                        let with_index = graph.get_lib_data(lib, *nid).unwrap();

                                        accumulator += with_index;
                                    }
                                }
                            }

                            let pos = self.tensors.len();
                            self.tensors.push(accumulator);

                            NetworkLeaf::LocalTensor(pos)
                        }
                    }
                    NetworkLeaf::LibraryKey(_) => {
                        let inds = graph.get_lib_data(lib, *nf).unwrap();
                        let mut accumulator = T::from(inds);
                        for (nid, t) in &targets[1..] {
                            match t {
                                NetworkLeaf::Scalar(_) => {
                                    return Err(TensorNetworkError::SumScalarTensor(
                                        "".to_string(),
                                    ));
                                }
                                NetworkLeaf::LocalTensor(t) => {
                                    accumulator += self.tensors[*t].refer();
                                }
                                NetworkLeaf::LibraryKey(_) => {
                                    let with = graph.get_lib_data(lib, *nid).unwrap();
                                    accumulator += with;
                                }
                            }
                        }

                        let pos = self.tensors.len();
                        self.tensors.push(accumulator);

                        NetworkLeaf::LocalTensor(pos)
                    }
                };

                graph.identify_nodes_without_self_edges(&all_nodes, NetworkNode::Leaf(new_node));
                Ok(graph)
            }
            NetworkOp::Function(f) => {
                let ops = graph
                    .graph
                    .iter_nodes()
                    .find(|(_, _, d)| matches!(d, NetworkNode::Op(NetworkOp::Function(_))));

                let (opid, children, _) = ops.unwrap();

                let mut child = None;
                for c in children {
                    if let Some(id) = graph.graph.involved_node_id(c)
                        && let NetworkNode::Leaf(l) = &graph.graph[id]
                    {
                        if let Some((nid, _)) = child {
                            if nid != id {
                                return Err(TensorNetworkError::Other(eyre!(
                                    "Cannot have more than one tensor argument to function"
                                )));
                            }
                        } else {
                            child = Some((id, l));
                        }
                    }
                }
                if let Some((child_id, leaf)) = child {
                    let new_node = match leaf {
                        NetworkLeaf::Scalar(s) => {
                            let s = self.scalar[*s].clone();
                            let pos = self.scalar.len();
                            let s = fn_lib.apply_scalar(&f, s)?;
                            self.scalar.push(s);

                            NetworkLeaf::Scalar(pos)
                        }
                        NetworkLeaf::LibraryKey(_) => {
                            let inds = graph.get_lib_data(lib, child_id).unwrap();
                            let t = fn_lib.apply(&f, T::from(inds))?;
                            let pos = self.tensors.len();
                            self.tensors.push(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        NetworkLeaf::LocalTensor(t) => {
                            let t = self.tensors[*t].clone();
                            let t = fn_lib.apply(&f, t)?;
                            let pos = self.tensors.len();
                            self.tensors.push(t);
                            NetworkLeaf::LocalTensor(pos)
                        }
                    };
                    graph.identify_nodes_without_self_edges(
                        &[child_id, opid],
                        NetworkNode::Leaf(new_node),
                    );
                    Ok(graph)
                } else {
                    Err(TensorNetworkError::ChildlessNeg)
                }
            }
            NetworkOp::Power(_) => {
                let mut pow = 0;
                let ops = graph.graph.iter_nodes().find(|(_, _, d)| {
                    if let NetworkNode::Op(NetworkOp::Power(i)) = d {
                        pow = *i;
                        true
                    } else {
                        false
                    }
                });

                let (opid, children, _) = ops.unwrap();

                let mut child = None;
                for c in children {
                    if let Some(id) = graph.graph.involved_node_id(c)
                        && let NetworkNode::Leaf(l) = &graph.graph[id]
                    {
                        if let Some((nid, _)) = child {
                            if nid != id {
                                return Err(TensorNetworkError::Other(eyre!(
                                    "Cannot have more than one tensor argument to power:{}",
                                    graph.dot()
                                )));
                            }
                        } else {
                            child = Some((id, l));
                        }
                    }
                }
                let n = pow.abs();
                if let Some((child_id, leaf)) = child {
                    let new_node = match leaf {
                        NetworkLeaf::Scalar(si) => {
                            if n == 0 {
                                NetworkLeaf::Scalar(*si)
                            } else {
                                let mut s = self.scalar[*si].clone();

                                for _ in 1..n {
                                    s *= self.scalar[*si].refer();
                                }

                                if pow < 0 {
                                    s = s.ref_one() / s;
                                }

                                let pos = self.scalar.len();
                                self.scalar.push(s);

                                NetworkLeaf::Scalar(pos)
                            }
                        }
                        NetworkLeaf::LibraryKey(a) => {
                            let inds = graph.get_lib_data(lib, child_id).unwrap();
                            let mut t = T::from(inds);

                            match pow {
                                0 => {
                                    let pos = self.scalar.len();
                                    let one = self.scalar[0].ref_one();
                                    self.scalar.push(one);
                                    NetworkLeaf::Scalar(pos)
                                }
                                1 => NetworkLeaf::LibraryKey(a.clone()),
                                _ => {
                                    let squares = n / 2;
                                    let mut square = t.contract(&t)?;

                                    if n % 2 == 1 {
                                        if n != 1 {
                                            for _ in 0..squares {
                                                square = square.contract(&square)?;
                                            }
                                            t = square.contract(&t)?;
                                        }

                                        if pow < 0 {
                                            if !t.is_scalar() {
                                                return Err(
                                                    TensorNetworkError::NegativeExponentNonScalar(
                                                        "".to_string(),
                                                    ),
                                                );
                                            } else {
                                                let mut s = Sc::from(t.scalar().unwrap());
                                                let pos = self.scalar.len();
                                                s = s.ref_one() / s;
                                                self.scalar.push(s);
                                                NetworkLeaf::Scalar(pos)
                                            }
                                        } else {
                                            let pos = self.tensors.len();
                                            self.tensors.push(t);
                                            NetworkLeaf::LocalTensor(pos)
                                        }
                                    } else {
                                        let mut s = Sc::from(square.scalar().unwrap());
                                        let sc = s.clone();
                                        for _ in 1..squares {
                                            s *= sc.refer();
                                        }
                                        let pos = self.scalar.len();
                                        if pow < 0 {
                                            s = s.ref_one() / s;
                                        }
                                        self.scalar.push(s);
                                        NetworkLeaf::Scalar(pos)
                                    }
                                }
                            }
                        }
                        NetworkLeaf::LocalTensor(ti) => {
                            let mut t = self.tensors[*ti].clone();
                            match pow {
                                0 => {
                                    let pos = self.scalar.len();
                                    let one = self.scalar[0].ref_one();
                                    self.scalar.push(one);
                                    NetworkLeaf::Scalar(pos)
                                }
                                1 => NetworkLeaf::LocalTensor(*ti),
                                _ => {
                                    let squares = n / 2;
                                    let mut square = t.contract(&t)?;

                                    if n % 2 == 1 {
                                        if n != 1 {
                                            for _ in 0..squares {
                                                square = square.contract(&square)?;
                                            }
                                            t = square.contract(&t)?;
                                        }
                                        if pow < 0 {
                                            if !t.is_scalar() {
                                                return Err(
                                                    TensorNetworkError::NegativeExponentNonScalar(
                                                        "".to_string(),
                                                    ),
                                                );
                                            } else {
                                                let mut s = Sc::from(t.scalar().unwrap());
                                                let pos = self.scalar.len();
                                                s = s.ref_one() / s;
                                                self.scalar.push(s);
                                                NetworkLeaf::Scalar(pos)
                                            }
                                        } else {
                                            let pos = self.tensors.len();
                                            self.tensors.push(t);
                                            NetworkLeaf::LocalTensor(pos)
                                        }
                                    } else {
                                        let mut s = Sc::from(square.scalar().unwrap());
                                        let sc = s.clone();
                                        for _ in 1..squares {
                                            s *= sc.refer();
                                        }
                                        let pos = self.scalar.len();
                                        if pow < 0 {
                                            s = s.ref_one() / s;
                                        }
                                        self.scalar.push(s);
                                        NetworkLeaf::Scalar(pos)
                                    }
                                }
                            }
                        }
                    };
                    graph.identify_nodes_without_self_edges(
                        &[child_id, opid],
                        NetworkNode::Leaf(new_node),
                    );
                    Ok(graph)
                } else {
                    Err(TensorNetworkError::ChildlessNeg)
                }
            }
        }
    }
}

pub trait Ref {
    type Ref<'a>
    where
        Self: 'a;
    fn refer(&self) -> Self::Ref<'_>;
}

impl Ref for f64 {
    type Ref<'a>
        = &'a f64
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

// #[cfg(feature = "shadowing")]
// pub mod levels;
#[cfg(feature = "shadowing")]
pub mod symbolica_interop;

#[cfg(test)]
mod tests;
