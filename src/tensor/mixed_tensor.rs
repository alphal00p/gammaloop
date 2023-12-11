use std::{collections::BTreeMap, ops::DerefMut};

use symbolica::{
    representations::{default::Linear, Atom, AtomBuilder, AtomSet},
    state::{BufferHandle, State, Workspace},
};

use crate::tensor::Slot;

use super::{
    ContractableWithDense, DenseTensor, HasTensorStructure, SparseTensor, VecSlotExtension,
};

pub struct SparseSymbolicTensor<'a, A: DerefMut<Target = Atom<P>>, P: AtomSet = Linear> {
    state: &'a State,
    workspace: &'a Workspace<P>,
    pub elements: BTreeMap<Vec<usize>, A>,
    pub structure: Vec<Slot>,
}

// impl<'a> ContractableWithDense<AtomBuilder<'a, BufferHandle<'a, Atom>>>
//     for SparseTensor<AtomBuilder<'a, BufferHandle<'a, Atom>>>
// {
//     fn contract_with_dense(
//         &self,
//         other: &DenseTensor<AtomBuilder<'a, BufferHandle<'a, Atom>>>,
//     ) -> Option<DenseTensor<AtomBuilder<'a, BufferHandle<'a, Atom>>>> {
//         if let Some((i, j)) = self.match_index(other) {
//             let final_structure = self.structure().merge_at(other.structure(), (i, j));
//             let mut result_data = vec![<$t>::default(); final_structure.size()];

//             let metric = self.structure()[i].representation.negative();

//             self.iter_fibers(i)
//                 .into_iter()
//                 .for_each(|(index_a, nonzeros, fiber_a)| {
//                     other
//                         .iter_fibers(j)
//                         .into_iter()
//                         .for_each(|(index_b, fiber_b)| {
//                             let result_index = final_structure
//                                 .flat_index(
//                                     &index_a[..i]
//                                         .iter()
//                                         .chain(&index_a[i + 1..])
//                                         .chain(&index_b[..j])
//                                         .chain(&index_b[j + 1..])
//                                         .cloned()
//                                         .collect::<Vec<_>>(),
//                                 )
//                                 .unwrap();
//                             for (i, k) in nonzeros.iter().enumerate() {
//                                 // Adjust indices for fetching from the other tensor
//                                 if metric[*k] {
//                                     result_data[result_index] -= (*fiber_a[i] * *fiber_b[*k]);
//                                 } else {
//                                     result_data[result_index] += *fiber_a[i] * *fiber_b[*k];
//                                 }
//                             }
//                         });
//                 });

//             let result = DenseTensor {
//                 data: result_data,
//                 structure: final_structure,
//             };

//             if result.traces().is_empty() {
//                 return Some(result);
//             } else {
//                 return Some(result.internal_contract());
//             }
//         }
//         None
//     }
// }
