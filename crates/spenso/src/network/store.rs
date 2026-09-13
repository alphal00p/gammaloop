use std::fmt::Debug;

use serde::{Deserialize, Serialize};

use crate::structure::slot::AbsInd;

use super::{
    graph::{NetworkGraph, ScalarRef},
    profile::{self, Counter, Timer},
};

pub trait TensorScalarStore: Default + TensorScalarStoreMapping {
    /// Tensor indices are local to the current graph and may change during execution.
    fn add_tensor(&mut self, tensor: Self::Tensor) -> usize;
    fn add_scalar(&mut self, scalar: Self::Scalar) -> usize;

    fn get_scalar(&self, index: usize) -> &Self::Scalar;
    fn get_scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        self.get_scalar(scalar.index())
    }
    fn get_tensor(&self, index: usize) -> &Self::Tensor;

    fn n_tensors(&self) -> usize;
    fn n_scalars(&self) -> usize;
    fn extend(&mut self, other: Self);
}

#[doc(hidden)]
pub trait NetworkStoreAccess {
    type Tensor;
    type Scalar;

    fn tensor(&self, index: usize) -> &Self::Tensor;
    fn scalar(&self, index: usize) -> &Self::Scalar;
    fn scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        self.scalar(scalar.index())
    }
    fn push_tensor(&mut self, tensor: Self::Tensor) -> usize;
    fn push_scalar(&mut self, scalar: Self::Scalar) -> usize;

    /// Drop unused tensors after installing a complete execution wave.
    /// The graph must include every pending reference to this store.
    fn retain_graph_tensors<K: Debug, FK: Debug, Aind: AbsInd>(
        &mut self,
        graph: &mut NetworkGraph<K, FK, Aind>,
    );
}

pub trait TensorScalarStoreMapping: Sized {
    type Tensor;

    type Scalar;
    type Store<T, S>: TensorScalarStoreMapping<Tensor = T, Scalar = S>;
    fn iter_scalars(&self) -> impl Iterator<Item = &Self::Scalar>;

    fn iter_tensors(&self) -> impl Iterator<Item = &Self::Tensor>;
    fn iter_scalars_mut(&mut self) -> impl Iterator<Item = &mut Self::Scalar>;
    fn iter_tensors_mut(&mut self) -> impl Iterator<Item = &mut Self::Tensor>;
    fn map<U, V>(
        self,
        scalar_map: impl FnMut(Self::Scalar) -> U,
        tensor_map: impl FnMut(Self::Tensor) -> V,
    ) -> Self::Store<V, U>;

    // fn map_self(
    //     self,
    //     scalar_map: impl FnMut(Self::Scalar) -> Self::Scalar,
    //     tensor_map: impl FnMut(Self::Tensor) -> Self::Tensor,
    // ) -> Self;

    fn map_result<U, V, Er>(
        self,
        scalar_map: impl FnMut(Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;

    // fn map_result_self<Er>(
    //     self,
    //     scalar_map: impl FnMut(Self::Scalar) -> Result<Self::Scalar, Er>,
    //     tensor_map: impl FnMut(Self::Tensor) -> Result<Self::Tensor, Er>,
    // ) -> Result<Self, Er>;

    fn map_ref<'a, U, V>(
        &'a self,
        scalar_map: impl FnMut(&'a Self::Scalar) -> U,
        tensor_map: impl FnMut(&'a Self::Tensor) -> V,
    ) -> Self::Store<V, U>;

    // fn map_ref_self(
    //     &self,
    //     scalar_map: impl FnMut(&Self::Scalar) -> Self::Scalar,
    //     tensor_map: impl FnMut(&Self::Tensor) -> Self::Tensor,
    // ) -> Self;

    fn map_ref_enumerate<U, V>(
        &self,
        scalar_map: impl FnMut((usize, &Self::Scalar)) -> U,
        tensor_map: impl FnMut((usize, &Self::Tensor)) -> V,
    ) -> Self::Store<V, U>;

    // fn map_ref_enumerate_self(
    //     &self,
    //     scalar_map: impl FnMut((usize, &Self::Scalar)) -> Self::Scalar,
    //     tensor_map: impl FnMut((usize, &Self::Tensor)) -> Self::Tensor,
    // ) -> Self;

    fn map_ref_result<U, V, Er>(
        &self,
        scalar_map: impl FnMut(&Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(&Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;

    // fn map_ref_result_self<Er>(
    //     &self,
    //     scalar_map: impl FnMut(&Self::Scalar) -> Result<Self::Scalar, Er>,
    //     tensor_map: impl FnMut(&Self::Tensor) -> Result<Self::Tensor, Er>,
    // ) -> Result<Self::Store<V, U>, Er>;

    fn map_ref_result_enumerate<U, V, Er>(
        &self,
        scalar_map: impl FnMut((usize, &Self::Scalar)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &Self::Tensor)) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;

    fn map_ref_mut<U, V>(
        &mut self,
        scalar_map: impl FnMut(&mut Self::Scalar) -> U,
        tensor_map: impl FnMut(&mut Self::Tensor) -> V,
    ) -> Self::Store<V, U>;

    fn map_ref_mut_enumerate<U, V>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut Self::Scalar)) -> U,
        tensor_map: impl FnMut((usize, &mut Self::Tensor)) -> V,
    ) -> Self::Store<V, U>;
    fn map_ref_mut_result<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut(&mut Self::Scalar) -> Result<U, Er>,
        tensor_map: impl FnMut(&mut Self::Tensor) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;

    fn map_ref_mut_result_enumerate<U, V, Er>(
        &mut self,
        scalar_map: impl FnMut((usize, &mut Self::Scalar)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &mut Self::Tensor)) -> Result<V, Er>,
    ) -> Result<Self::Store<V, U>, Er>;
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
)]
pub struct NetworkStore<T, S> {
    pub tensors: Vec<T>,
    // pub params: AHashSet<Atom>,
    pub scalar: Vec<S>,
    pub scalar_aliases: Vec<Option<S>>,
}

#[doc(hidden)]
pub struct NetworkStoreOverlay<'a, T, S> {
    base_tensors: &'a [T],
    base_scalars: &'a [S],
    base_scalar_aliases: &'a [Option<S>],
    pub tensors: Vec<T>,
    pub scalar: Vec<S>,
    pub scalar_aliases: Vec<Option<S>>,
}

impl<'a, T, S> NetworkStoreOverlay<'a, T, S> {
    pub fn new(base: &'a NetworkStore<T, S>) -> Self {
        Self {
            base_tensors: &base.tensors,
            base_scalars: &base.scalar,
            base_scalar_aliases: &base.scalar_aliases,
            tensors: Vec::new(),
            scalar: Vec::new(),
            scalar_aliases: Vec::new(),
        }
    }

    pub fn into_additions(self) -> (Vec<T>, Vec<S>, Vec<Option<S>>) {
        (self.tensors, self.scalar, self.scalar_aliases)
    }
}

impl<T, S> Default for NetworkStore<T, S> {
    fn default() -> Self {
        NetworkStore {
            tensors: vec![],
            scalar: vec![],
            scalar_aliases: vec![],
        }
    }
}

impl<T, S> NetworkStoreAccess for NetworkStore<T, S> {
    type Tensor = T;
    type Scalar = S;

    fn tensor(&self, index: usize) -> &Self::Tensor {
        &self.tensors[index]
    }

    fn scalar(&self, index: usize) -> &Self::Scalar {
        &self.scalar[index]
    }

    fn scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        match scalar {
            ScalarRef::Store(index) => &self.scalar[index],
            ScalarRef::Alias(index) => self
                .scalar_aliases
                .get(index)
                .and_then(Option::as_ref)
                .unwrap_or(&self.scalar[index]),
        }
    }

    fn push_tensor(&mut self, tensor: Self::Tensor) -> usize {
        let index = self.tensors.len();
        self.tensors.push(tensor);
        index
    }

    fn push_scalar(&mut self, scalar: Self::Scalar) -> usize {
        let index = self.scalar.len();
        self.scalar.push(scalar);
        self.scalar_aliases.push(None);
        index
    }

    fn retain_graph_tensors<K: Debug, FK: Debug, Aind: AbsInd>(
        &mut self,
        graph: &mut NetworkGraph<K, FK, Aind>,
    ) {
        if self.tensors.is_empty() {
            return;
        }
        // MAX marks an unreferenced slot, never a graph index. Mark every
        // reference, including shared IDs, before moving any tensor.
        let mut indices = vec![usize::MAX; self.tensors.len()];
        graph.map_tensor_refs(|index| {
            indices[index] = index;
            index
        });
        let mut retained = 0;
        let mut old_indices = indices.iter_mut();
        self.tensors.retain(|_| {
            let index = old_indices.next().unwrap();
            if *index == usize::MAX {
                return false;
            }
            *index = retained;
            retained += 1;
            true
        });
        graph.map_tensor_refs(|index| indices[index]);
    }
}

impl<T, S> NetworkStoreAccess for NetworkStoreOverlay<'_, T, S> {
    type Tensor = T;
    type Scalar = S;

    fn tensor(&self, index: usize) -> &Self::Tensor {
        if index < self.base_tensors.len() {
            &self.base_tensors[index]
        } else {
            &self.tensors[index - self.base_tensors.len()]
        }
    }

    fn scalar(&self, index: usize) -> &Self::Scalar {
        if index < self.base_scalars.len() {
            &self.base_scalars[index]
        } else {
            &self.scalar[index - self.base_scalars.len()]
        }
    }

    fn scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        let index = scalar.index();
        match scalar {
            ScalarRef::Store(_) => self.scalar(index),
            ScalarRef::Alias(_) => {
                if index < self.base_scalars.len() {
                    self.base_scalar_aliases
                        .get(index)
                        .and_then(Option::as_ref)
                        .unwrap_or(&self.base_scalars[index])
                } else {
                    let local = index - self.base_scalars.len();
                    self.scalar_aliases
                        .get(local)
                        .and_then(Option::as_ref)
                        .unwrap_or(&self.scalar[local])
                }
            }
        }
    }

    fn push_tensor(&mut self, tensor: Self::Tensor) -> usize {
        let index = self.base_tensors.len() + self.tensors.len();
        self.tensors.push(tensor);
        index
    }

    fn push_scalar(&mut self, scalar: Self::Scalar) -> usize {
        let index = self.base_scalars.len() + self.scalar.len();
        self.scalar.push(scalar);
        self.scalar_aliases.push(None);
        index
    }

    fn retain_graph_tensors<K: Debug, FK: Debug, Aind: AbsInd>(
        &mut self,
        _graph: &mut NetworkGraph<K, FK, Aind>,
    ) {
        // Workers borrow the shared base and can hold references outside their
        // operation graph. Reclamation belongs to the joined, merged wave.
    }
}

impl<T, S> TensorScalarStore for NetworkStore<T, S> {
    fn n_tensors(&self) -> usize {
        self.tensors.len()
    }

    fn n_scalars(&self) -> usize {
        self.scalar.len()
    }

    fn extend(&mut self, other: Self) {
        let _span = profile::span(Timer::StoreExtend);
        profile::bump(Counter::StoreExtend, 1);
        self.tensors.extend(other.tensors);
        self.scalar.extend(other.scalar);
        self.scalar_aliases.extend(other.scalar_aliases);
    }

    fn add_scalar(&mut self, scalar: Self::Scalar) -> usize {
        let id = self.scalar.len();
        self.scalar.push(scalar);
        self.scalar_aliases.push(None);
        id
    }

    fn add_tensor(&mut self, tensor: Self::Tensor) -> usize {
        let id = self.tensors.len();
        self.tensors.push(tensor);
        id
    }

    fn get_scalar(&self, index: usize) -> &Self::Scalar {
        &self.scalar[index]
    }

    fn get_scalar_ref(&self, scalar: ScalarRef) -> &Self::Scalar {
        match scalar {
            ScalarRef::Store(index) => &self.scalar[index],
            ScalarRef::Alias(index) => self
                .scalar_aliases
                .get(index)
                .and_then(Option::as_ref)
                .unwrap_or(&self.scalar[index]),
        }
    }

    fn get_tensor(&self, index: usize) -> &Self::Tensor {
        &self.tensors[index]
    }
}

impl<T, S> TensorScalarStoreMapping for NetworkStore<T, S> {
    type Store<U, V> = NetworkStore<U, V>;
    type Scalar = S;
    type Tensor = T;

    fn iter_scalars(&self) -> impl Iterator<Item = &Self::Scalar> {
        self.scalar.iter()
    }

    fn iter_tensors(&self) -> impl Iterator<Item = &Self::Tensor> {
        self.tensors.iter()
    }

    fn iter_scalars_mut(&mut self) -> impl Iterator<Item = &mut Self::Scalar> {
        self.scalar.iter_mut()
    }

    fn iter_tensors_mut(&mut self) -> impl Iterator<Item = &mut Self::Tensor> {
        self.tensors.iter_mut()
    }

    fn map<U, V>(
        self,
        mut scalar_map: impl FnMut(S) -> U,
        tensor_map: impl FnMut(T) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.into_iter().map(tensor_map).collect(),
            scalar: self.scalar.into_iter().map(&mut scalar_map).collect(),
            scalar_aliases: self
                .scalar_aliases
                .into_iter()
                .map(|alias| alias.map(&mut scalar_map))
                .collect(),
        }
    }

    fn map_result<U, V, Er>(
        self,
        mut scalar_map: impl FnMut(S) -> Result<U, Er>,
        tensor_map: impl FnMut(T) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .into_iter()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .into_iter()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .into_iter()
                .map(|alias| alias.map(&mut scalar_map).transpose())
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref<'a, U, V>(
        &'a self,
        mut scalar_map: impl FnMut(&'a S) -> U,
        tensor_map: impl FnMut(&'a T) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.iter().map(tensor_map).collect(),
            scalar: self.scalar.iter().map(&mut scalar_map).collect(),
            scalar_aliases: self
                .scalar_aliases
                .iter()
                .map(|alias| alias.as_ref().map(&mut scalar_map))
                .collect(),
        }
    }

    fn map_ref_enumerate<U, V>(
        &self,
        mut scalar_map: impl FnMut((usize, &S)) -> U,
        tensor_map: impl FnMut((usize, &T)) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.iter().enumerate().map(tensor_map).collect(),
            scalar: self
                .scalar
                .iter()
                .enumerate()
                .map(&mut scalar_map)
                .collect(),
            scalar_aliases: self
                .scalar_aliases
                .iter()
                .enumerate()
                .map(|(index, alias)| alias.as_ref().map(|scalar| scalar_map((index, scalar))))
                .collect(),
        }
    }

    fn map_ref_result<U, V, Er>(
        &self,
        mut scalar_map: impl FnMut(&S) -> Result<U, Er>,
        tensor_map: impl FnMut(&T) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .iter()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .iter()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .iter()
                .map(|alias| alias.as_ref().map(&mut scalar_map).transpose())
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref_result_enumerate<U, V, Er>(
        &self,
        mut scalar_map: impl FnMut((usize, &S)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &T)) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .iter()
                .enumerate()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .iter()
                .enumerate()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .iter()
                .enumerate()
                .map(|(index, alias)| {
                    alias
                        .as_ref()
                        .map(|scalar| scalar_map((index, scalar)))
                        .transpose()
                })
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref_mut<U, V>(
        &mut self,
        mut scalar_map: impl FnMut(&mut S) -> U,
        tensor_map: impl FnMut(&mut T) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self.tensors.iter_mut().map(tensor_map).collect(),
            scalar: self.scalar.iter_mut().map(&mut scalar_map).collect(),
            scalar_aliases: self
                .scalar_aliases
                .iter_mut()
                .map(|alias| alias.as_mut().map(&mut scalar_map))
                .collect(),
        }
    }

    fn map_ref_mut_enumerate<U, V>(
        &mut self,
        mut scalar_map: impl FnMut((usize, &mut S)) -> U,
        tensor_map: impl FnMut((usize, &mut T)) -> V,
    ) -> NetworkStore<V, U> {
        NetworkStore {
            tensors: self
                .tensors
                .iter_mut()
                .enumerate()
                .map(tensor_map)
                .collect(),
            scalar: self
                .scalar
                .iter_mut()
                .enumerate()
                .map(&mut scalar_map)
                .collect(),
            scalar_aliases: self
                .scalar_aliases
                .iter_mut()
                .enumerate()
                .map(|(index, alias)| alias.as_mut().map(|scalar| scalar_map((index, scalar))))
                .collect(),
        }
    }

    fn map_ref_mut_result<U, V, Er>(
        &mut self,
        mut scalar_map: impl FnMut(&mut S) -> Result<U, Er>,
        tensor_map: impl FnMut(&mut T) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .iter_mut()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .iter_mut()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .iter_mut()
                .map(|alias| alias.as_mut().map(&mut scalar_map).transpose())
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }

    fn map_ref_mut_result_enumerate<U, V, Er>(
        &mut self,
        mut scalar_map: impl FnMut((usize, &mut S)) -> Result<U, Er>,
        tensor_map: impl FnMut((usize, &mut T)) -> Result<V, Er>,
    ) -> Result<NetworkStore<V, U>, Er> {
        Ok(NetworkStore {
            tensors: self
                .tensors
                .iter_mut()
                .enumerate()
                .map(tensor_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar: self
                .scalar
                .iter_mut()
                .enumerate()
                .map(&mut scalar_map)
                .collect::<Result<Vec<_>, Er>>()?,
            scalar_aliases: self
                .scalar_aliases
                .iter_mut()
                .enumerate()
                .map(|(index, alias)| {
                    alias
                        .as_mut()
                        .map(|scalar| scalar_map((index, scalar)))
                        .transpose()
                })
                .collect::<Result<Vec<_>, Er>>()?,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Mutex;

    use linnet::{
        half_edge::{builder::HedgeGraphBuilder, involution::Flow},
        permutation::Permutation,
    };

    use crate::{
        network::graph::{NetworkEdge, NetworkLeaf, NetworkNode, ScaledTensorRef},
        structure::PermutedStructure,
    };

    use super::{
        NetworkGraph, NetworkStore, NetworkStoreAccess, NetworkStoreOverlay, ScalarRef,
        TensorScalarStore,
    };

    #[test]
    fn tensor_retention_moves_shared_payloads_and_preserves_every_leaf_variant() {
        let library = NetworkLeaf::LibraryKey {
            key: PermutedStructure {
                structure: 7i8,
                rep_permutation: Permutation::id(0),
                index_permutation: Permutation::id(0),
            },
            indices: vec![],
        };
        let mut builder = HedgeGraphBuilder::new();
        for leaf in [
            NetworkLeaf::LocalTensor(2),
            NetworkLeaf::TensorSum(vec![2, 5, 2]),
            NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled_ref(8, ScalarRef::Alias(0))),
            NetworkLeaf::ScaledTensorSum(vec![
                ScaledTensorRef::tensor(5),
                ScaledTensorRef::scaled(8, 1),
            ]),
            NetworkLeaf::Scalar(ScalarRef::Alias(0)),
            library.clone(),
        ] {
            let node = builder.add_node(NetworkNode::Leaf(leaf));
            builder.add_external_edge(node, NetworkEdge::Head, true, Flow::Source);
        }
        let mut graph: NetworkGraph<i8, i8> = builder.into();
        // Mutex deliberately has no Clone implementation: compaction must move payloads.
        let mut store = NetworkStore {
            tensors: (0..10).map(Mutex::new).collect(),
            scalar: vec![11, 13],
            scalar_aliases: vec![Some(17), None],
        };
        let expected = vec![
            NetworkLeaf::LocalTensor(0),
            NetworkLeaf::TensorSum(vec![0, 1, 0]),
            NetworkLeaf::ScaledTensor(ScaledTensorRef::scaled_ref(2, ScalarRef::Alias(0))),
            NetworkLeaf::ScaledTensorSum(vec![
                ScaledTensorRef::tensor(1),
                ScaledTensorRef::scaled(2, 1),
            ]),
            NetworkLeaf::Scalar(ScalarRef::Alias(0)),
            library,
        ];
        for _ in 0..2 {
            store.retain_graph_tensors(&mut graph);
            assert_eq!(
                store
                    .tensors
                    .iter()
                    .map(|tensor| *tensor.lock().unwrap())
                    .collect::<Vec<_>>(),
                [2, 5, 8]
            );
            assert_eq!(
                graph
                    .graph
                    .iter_nodes()
                    .map(|(_, _, node)| {
                        let NetworkNode::Leaf(leaf) = node else {
                            panic!("expected leaf")
                        };
                        leaf.clone()
                    })
                    .collect::<Vec<_>>(),
                expected
            );
            assert_eq!(store.scalar, [11, 13]);
            assert_eq!(store.scalar_aliases, [Some(17), None]);
        }
        let mut scalar_graph: NetworkGraph<i8, i8> = NetworkGraph::scalar(0);
        store.retain_graph_tensors(&mut scalar_graph);
        assert!(store.tensors.is_empty());
        assert_eq!(store.scalar_ref(ScalarRef::Alias(0)), &17);
    }

    #[test]
    fn overlay_retention_keeps_base_and_unmerged_addition_indices() {
        let base = NetworkStore {
            tensors: vec![10, 20, 30],
            scalar: vec![7],
            scalar_aliases: vec![Some(11)],
        };
        let mut overlay = NetworkStoreOverlay::new(&base);
        assert_eq!(overlay.push_tensor(40), 3);
        assert_eq!(overlay.push_tensor(50), 4);
        let mut graph: NetworkGraph<i8, i8> = NetworkGraph::scalar(0);
        let root = graph.graph.node_id(graph.head());
        graph.graph[root] = NetworkNode::Leaf(NetworkLeaf::TensorSum(vec![1, 4]));
        let before = graph.clone();
        overlay.retain_graph_tensors(&mut graph);
        assert_eq!(graph, before);
        assert_eq!(overlay.tensor(1), &20);
        assert_eq!(overlay.tensor(4), &50);
        assert_eq!(overlay.scalar_ref(ScalarRef::Alias(0)), &11);
        assert_eq!(overlay.into_additions(), (vec![40, 50], vec![], vec![]));
        assert_eq!(base.tensors, [10, 20, 30]);
    }

    #[cfg(feature = "shadowing")]
    #[test]
    fn tensor_retention_preserves_unreferenced_alias_definitions() {
        use symbolica::{atom::Atom, parse};

        use crate::network::{Network, tags::scalar_store_alias};

        let mut network: Network<NetworkStore<Mutex<i32>, Atom>, i8, i8> =
            Network::from_scalar(parse!("x+y"));
        let nested = network
            .store
            .add_scalar(scalar_store_alias(0) + parse!("z"));
        let aliases = network.alias_scalar_refs(|_, _| true);
        network.store.add_tensor(Mutex::new(1));
        network.store.retain_graph_tensors(&mut network.graph);
        assert!(network.store.tensors.is_empty());
        assert_eq!(
            network.resolve_scalar_aliases(&aliases, scalar_store_alias(nested)),
            parse!("x+y+z")
        );
        assert_eq!(aliases.aliased_indices().collect::<Vec<_>>(), [0, 1]);
    }

    #[test]
    fn completed_execution_waves_reclaim_tensor_inputs_in_every_strategy() {
        use crate::{
            network::{
                ExecutionResult, Network, Sequential, SequentialExtract, SequentialRef,
                SmallestDegree,
                library::{DummyKey, DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
            },
            structure::{
                OrderedStructure,
                representation::{Euclidean, RepName},
            },
            tensors::data::DenseTensor,
        };

        #[cfg(feature = "shadowing")]
        let _parallel = crate::symbolic_parallelism::scoped_symbolica_rayon_setting_for_test(
            crate::symbolic_parallelism::SymbolicParallelism::Parallel,
            || true,
        );
        type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
        type Net = Network<NetworkStore<Tensor, f64>, DummyKey, DummyKey>;
        type LibTensor = DummyLibraryTensor<Tensor>;
        type Lib = DummyLibrary<Tensor, DummyKey>;
        type FnLib = ErroringLibrary<DummyKey>;
        let lib = Lib::new();
        let fn_lib = FnLib::new();
        let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).structure;
        let tensor = |data| Net::from_tensor(Tensor::from_data(data, structure.clone()).unwrap());
        // The two inner sums can run in separate overlays before later waves
        // negate and combine them. Every stage must use the remapped indices.
        let original = -(tensor(vec![1.0, 2.0]) + tensor(vec![3.0, 4.0]))
            + -(tensor(vec![5.0, 6.0]) + tensor(vec![7.0, 8.0]));
        assert_eq!(original.store.tensors.len(), 4);
        for strategy in 0..4 {
            let mut network = original.clone();
            match strategy {
                0 => network
                    .execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib),
                1 => network
                    .execute::<SequentialRef, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib),
                2 => network.execute::<SequentialExtract, SmallestDegree, LibTensor, Lib, FnLib>(
                    &lib, &fn_lib,
                ),
                _ => {
                    network.execute_parallel::<SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
                }
            }
            .unwrap();
            assert_eq!(network.store.tensors.len(), 1);
            let (NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)), _, _) =
                network.graph.result().unwrap()
            else {
                panic!("expected a materialized tensor")
            };
            assert_eq!(*index, 0);
            assert_eq!(network.store.tensors[0].data, [-16.0, -20.0]);
        }

        // A self-loop-only network finishes before scheduling an operation wave.
        let structure = OrderedStructure::new(vec![
            Euclidean {}.new_slot(2, 1),
            Euclidean {}.new_slot(2, 1),
        ])
        .structure;
        let mut trace =
            Net::from_tensor(Tensor::from_data(vec![1.0, 2.0, 3.0, 4.0], structure).unwrap());
        trace
            .execute::<Sequential, SmallestDegree, LibTensor, Lib, FnLib>(&lib, &fn_lib)
            .unwrap();
        assert_eq!(trace.store.tensors.len(), 1);
        let ExecutionResult::Val(value) = trace.result_scalar().unwrap() else {
            panic!("expected a trace")
        };
        assert_eq!(*value, 5.0);
    }
}
