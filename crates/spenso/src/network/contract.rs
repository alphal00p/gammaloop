use std::{
    fmt::{Debug, Display},
    ops::MulAssign,
};

use crate::{
    algebra::ScalarMul,
    contraction::Contract,
    network::graph::{NetworkLeaf, NetworkNode, NetworkOp},
    structure::{
        HasStructure, PermutedStructure, TensorStructure,
        permuted::PermuteTensor,
        slot::{AbsInd, IsAbstractSlot},
    },
};

use super::{
    Ref, TensorNetworkError,
    graph::NetworkGraph,
    library::{Library, LibraryTensor},
    store::NetworkStore,
};

pub struct SmallestDegree;

pub struct SmallestDegreeIter<const N: usize>;

pub struct ContractScalars;

pub struct SingleSmallestDegree<const D: bool>;
pub trait ContractionStrategy<E, L, K, FK, Aind>: Sized {
    #[allow(clippy::result_large_err, clippy::type_complexity)]
    fn contract(
        executor: &mut E,
        graph: NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<(NetworkGraph<K, FK, Aind>, bool), TensorNetworkError<K, FK>>
    where
        K: Display,
        FK: Display;
}

impl<
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
> ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind> for ContractScalars
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    /// Contract all scalars into One,
    /// If there is a *single* other tensor, then scalar multiply with this tensor and then remove the op node
    /// If there are more than one other tensors,then only reduce to a single scalar, keeping the op node
    /// If there are no other tensors, then reduce to a single scalar and remove the op node
    fn contract(
        executor: &mut NetworkStore<T, Sc>,
        mut graph: NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<(NetworkGraph<K, FK, Aind>, bool), TensorNetworkError<K, FK>>
    where
        K: Display,
    {
        graph.sync_order();
        let mut other = None;
        let mut remove_op_node = true;
        let mut head = None;
        // println!("{}", graph.dot());

        let (mut scalars, mut scalar_nodes): (Vec<_>, Vec<_>) = graph
            .graph
            .iter_nodes()
            .filter_map(|(nid, _, c)| {
                if let NetworkNode::Leaf(l) = c {
                    match l {
                        NetworkLeaf::Scalar(i) => Some((*i, nid)),
                        _ => {
                            if other.is_none() {
                                other = Some(nid);
                            } else {
                                remove_op_node = false;
                            }
                            None
                        }
                    }
                } else {
                    if let NetworkNode::Op(NetworkOp::Product) = c {
                        if head.is_some() {
                            panic!("multiple heads")
                        }
                        head = Some(nid);
                    }
                    None
                }
            })
            .collect();

        // println!("Scalars {scalars:?} nodes {scalar_nodes:?}");

        if let Some(f) = scalars.pop() {
            let mut acc = executor.scalar[f].clone();
            let mut more_than_one_scalar = false;

            // Accumulate all scalars
            for si in scalars {
                more_than_one_scalar = true;
                acc *= executor.scalar[si].refer();
            }

            let new_node = if remove_op_node {
                if let Some(head) = head {
                    //Should always have an op node
                    // Here we add it to the list of nodes we want to merge
                    scalar_nodes.push(head);
                }
                if let Some(other) = other {
                    // Since remove_op_node is true,
                    // Then this means other is the single other tensor node
                    // We now perform scalar multiplication, creating a new tensor node
                    //
                    // println!("Single  tensor {other}, scalar multiplying");

                    scalar_nodes.push(other);
                    let NetworkNode::Leaf(l) = &graph.graph[other] else {
                        unreachable!("aa")
                    };
                    match l {
                        NetworkLeaf::LocalTensor(l) => {
                            let a = executor.tensors[*l].scalar_mul(&acc).unwrap();
                            let pos = executor.tensors.len();
                            executor.tensors.push(a);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        NetworkLeaf::LibraryKey(_) => {
                            let inds = graph.get_lib_data(lib, other).unwrap();
                            let a = inds.scalar_mul(&acc).unwrap();

                            let pos = executor.tensors.len();
                            executor.tensors.push(a);
                            NetworkLeaf::LocalTensor(pos)
                        }
                        _ => {
                            unreachable!("aa")
                        }
                    }
                } else {
                    // This means that we only have scalars,
                    // We create the resulting scalar node, and return that
                    // println!("Only scalars");
                    let pos = executor.scalar.len();
                    executor.scalar.push(acc);
                    NetworkLeaf::Scalar(pos)
                }
            } else {
                // This means that we have multiple tensors, we thus will only merge the scalars
                // We create the resulting scalar node, and return that
                // println!("Multiple tensors only merge scalars");
                if more_than_one_scalar {
                    let pos = executor.scalar.len();
                    executor.scalar.push(acc);
                    NetworkLeaf::Scalar(pos)
                } else {
                    // println!("Single scalar");
                    return Ok((graph, false));
                }
            };

            if !remove_op_node {
                graph.identify_nodes_without_self_edges_merge_heads(
                    &scalar_nodes,
                    NetworkNode::Leaf(new_node),
                );
            } else {
                graph.identify_nodes_without_self_edges(&scalar_nodes, NetworkNode::Leaf(new_node));
            }
            // println!("{}", graph.dot());
            Ok((graph, true))
        } else {
            let mut didsmth = false;
            if remove_op_node
                && let Some(other) = other
                && let Some(head) = head
            {
                let v = graph.graph[other].clone();
                graph.identify_nodes_without_self_edges(&[head, other], v);
                didsmth = true;
            }

            Ok((graph, didsmth))
        }
    }
}

impl<
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
> ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind> for SmallestDegree
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut NetworkStore<T, Sc>,
        graph: NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<(NetworkGraph<K, FK, Aind>, bool), TensorNetworkError<K, FK>>
    where
        K: Display,
    {
        // println!("Contracting scalars");
        let (mut graph, mut didsmth) = ContractScalars::contract(executor, graph, lib)?;

        // println!("Contracted scalars");

        while {
            let (newgraph, smth) = SingleSmallestDegree::<false>::contract(executor, graph, lib)?;
            graph = newgraph;
            smth
        } {
            didsmth |= true
        }

        let (graph, _) = ContractScalars::contract(executor, graph, lib)?;

        Ok((graph, didsmth))
    }
}

impl<
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
    const N: usize,
> ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind> for SmallestDegreeIter<N>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut NetworkStore<T, Sc>,
        graph: NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<(NetworkGraph<K, FK, Aind>, bool), TensorNetworkError<K, FK>>
    where
        K: Display,
    {
        let (mut graph, mut didsmth) = ContractScalars::contract(executor, graph, lib)?;

        for _ in 0..N {
            let (newgraph, smth) = SingleSmallestDegree::<false>::contract(executor, graph, lib)?;
            graph = newgraph;
            didsmth |= smth;
        }

        let (graph, _) = ContractScalars::contract(executor, graph, lib)?;

        Ok((graph, didsmth))
    }
}

impl<
    LT: LibraryTensor + Clone,
    T: HasStructure
        + TensorStructure
        + Clone
        + Contract<LCM = T>
        + ScalarMul<Sc, Output = T>
        + Contract<LT::WithIndices, LCM = T>
        + From<LT::WithIndices>,
    L: Library<T::Structure, Key = K, Value = PermutedStructure<LT>>,
    Sc: for<'a> MulAssign<Sc::Ref<'a>>
        + Clone
        + for<'a> MulAssign<T::ScalarRef<'a>>
        + From<T::Scalar>
        + Ref,
    K: Display + Debug + Clone,
    FK: Display + Debug + Clone,
    Aind: AbsInd,
    const D: bool,
> ContractionStrategy<NetworkStore<T, Sc>, L, K, FK, Aind> for SingleSmallestDegree<D>
where
    LT::WithIndices: Contract<LT::WithIndices, LCM = T>
        + ScalarMul<Sc, Output = T>
        + PermuteTensor<Permuted = LT::WithIndices>,
    <LT::WithIndices as HasStructure>::Structure: Display,
    T::Structure: Display,
    <<LT::WithIndices as HasStructure>::Structure as TensorStructure>::Slot:
        IsAbstractSlot<Aind = Aind>,
{
    fn contract(
        executor: &mut NetworkStore<T, Sc>,
        mut graph: NetworkGraph<K, FK, Aind>,
        lib: &L,
    ) -> Result<(NetworkGraph<K, FK, Aind>, bool), TensorNetworkError<K, FK>>
    where
        K: Display,
    {
        graph.sync_order();
        if D {
            println!("Contracting {}", graph.dot());
        }

        let mut last_tensor = None;
        let edge_to_contract = graph
            .graph
            .iter_nodes()
            .filter(|(_, _, d)| d.is_tensor())
            .filter_map(|(nid1, a, n1)| {
                let mut degree = 0;
                let mut first = None;
                for h in a {
                    if graph.graph[[&h]].is_slot() && graph.graph.inv(h) != h {
                        first = Some(h); //only contract slot hedges
                        degree += 1
                    }
                }

                let nid2 = if degree == 0 {
                    //no internal slots to contract
                    // contract with last tensor (give max  weight  so only happens when there are no internal slots)
                    degree = i32::MAX;
                    if let Some(last_tensor) = last_tensor {
                        last_tensor
                    } else {
                        last_tensor = Some(nid1);
                        return None;
                    }
                } else {
                    graph.graph.involved_node_id(first?)?
                };

                let n2 = &graph.graph[nid2];

                last_tensor = Some(nid1);

                Some((degree, nid1, n1, nid2, n2))
            })
            .min_by_key(|(degree, _, _, _, _)| *degree);

        if let Some((_, nid1, n1, nid2, n2)) = edge_to_contract {
            if D {
                println!("Contracting {} with {}", nid1, nid2);
            }
            let new_node = match (n1, n2) {
                (NetworkNode::Leaf(_), NetworkNode::Op(NetworkOp::Product))
                | (NetworkNode::Op(NetworkOp::Product), NetworkNode::Leaf(_)) => {
                    return Err(TensorNetworkError::SlotEdgeToProdNode);
                }
                (NetworkNode::Leaf(l1), NetworkNode::Leaf(l2)) => match (l1, l2) {
                    (NetworkLeaf::Scalar(_), _) | (_, NetworkLeaf::Scalar(_)) => {
                        return Err(TensorNetworkError::SlotEdgeToScalarNode);
                    }

                    (NetworkLeaf::LocalTensor(l1), NetworkLeaf::LocalTensor(l2)) => {
                        if D {
                            let st1 = executor.tensors[*l1].structure();
                            let st2 = executor.tensors[*l2].structure();

                            println!("Contracting {} with {}", st1, st2);
                        }

                        let contracted = executor.tensors[*l1].contract(&executor.tensors[*l2])?;

                        if D {
                            println!("Obtained {}", contracted.structure());
                        }
                        let pos = executor.tensors.len();
                        executor.tensors.push(contracted);

                        NetworkLeaf::LocalTensor(pos)
                    }
                    (NetworkLeaf::LibraryKey(_), NetworkLeaf::LocalTensor(l2)) => {
                        let l1 = graph.get_lib_data(lib, nid1).unwrap();
                        if D {
                            let st1 = l1.structure();
                            let st2 = executor.tensors[*l2].structure();
                            println!("Contracting {} with {}", st1, st2);
                        }

                        let contracted = executor.tensors[*l2].contract(&l1)?;
                        if D {
                            println!("Obtained {}", contracted.structure());
                        }
                        let pos = executor.tensors.len();
                        executor.tensors.push(contracted);
                        NetworkLeaf::LocalTensor(pos)
                    }

                    (NetworkLeaf::LocalTensor(l2), NetworkLeaf::LibraryKey(_)) => {
                        let l1 = graph.get_lib_data(lib, nid2).unwrap();
                        if D {
                            let st1 = l1.structure();
                            let st2 = executor.tensors[*l2].structure();
                            println!("Contracting {} with {}", st2, st1);
                        }

                        let contracted = executor.tensors[*l2].contract(&l1)?;
                        if D {
                            println!("Obtained {}", contracted.structure());
                        }
                        let pos = executor.tensors.len();
                        executor.tensors.push(contracted);

                        NetworkLeaf::LocalTensor(pos)
                    }
                    (NetworkLeaf::LibraryKey(_), NetworkLeaf::LibraryKey(_)) => {
                        let l1 = graph.get_lib_data(lib, nid1).unwrap();

                        let l2 = graph.get_lib_data(lib, nid2).unwrap();
                        if D {
                            let st1 = l1.structure();
                            let st2 = l2.structure();
                            println!("Contracting {} with {}", st2, st1);
                        }

                        let contracted = l1.contract(&l2)?;
                        if D {
                            println!("Obtained {}", contracted.structure());
                        }
                        let pos = executor.tensors.len();
                        executor.tensors.push(contracted);

                        NetworkLeaf::LocalTensor(pos)
                    }
                },
                (a, b) => {
                    return Err(TensorNetworkError::CannotContractEdgeBetween(
                        a.clone(),
                        b.clone(),
                    ));
                }
            };
            graph.identify_nodes_without_self_edges_merge_heads(
                &[nid1, nid2],
                NetworkNode::Leaf(new_node),
            );
            Ok((graph, true))
        } else {
            Ok((graph, false))
        }
    }
}
