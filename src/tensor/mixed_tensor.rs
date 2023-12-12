use std::{collections::BTreeMap, ops::DerefMut};

use symbolica::{
    representations::{default::Linear, AsAtomView, Atom, AtomBuilder, AtomSet},
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
impl SparseTensor<Atom> {
    fn contract_with_dense(
        &self,
        other: &DenseTensor<Atom>,
        workspace: &Workspace,
        state: &State,
    ) -> Option<DenseTensor<Atom>> {
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let zero = workspace.new_num(0);

            let neutral_summand = zero.builder(&state, &workspace);
            let mut result_data = (0..final_structure.size())
                .map(|_| zero.builder(&state, &workspace))
                .collect::<Vec<_>>();
            let metric = self.structure()[i].representation.negative();

            for (index_a, nonzeros, fiber_a) in self.iter_fibers(i) {
                for (index_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = final_structure
                        .flat_index(
                            &index_a[..i]
                                .iter()
                                .chain(&index_a[i + 1..])
                                .chain(&index_b[..j])
                                .chain(&index_b[j + 1..])
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                        .unwrap();
                    for (i, k) in nonzeros.iter().enumerate() {
                        if metric[*k] {
                            result_data[result_index] = -(fiber_a[i].builder(state, workspace)
                                * fiber_b[*k])
                                + result_data[result_index].as_atom_view();
                        }
                    }
                }
            }
        }
        None
    }
}
