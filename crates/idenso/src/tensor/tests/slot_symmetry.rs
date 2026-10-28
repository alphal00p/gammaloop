use crate::tensor::SymbolicTensor;
use crate::tensor::composition::MatrixChannel;
use spenso::network::tensor_symmetry::TensorSlotSymmetry;
use spenso::shadowing::FactorProjector;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
use spenso::structure::slot::IsAbstractSlot;
use spenso::tensor_symbol;
use std::collections::HashMap;
use symbolica::atom::{Atom, AtomCore, NamespacedSymbol, Symbol, SymbolAttribute, UserData};
use symbolica::function;

type Tensor = SymbolicTensor<PartialStructure>;
fn slot(n: usize) -> Atom {
    LibraryRep::from(Minkowski {})
        .new_rep(4)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(n))
        .to_atom()
}
fn head(name: &str) -> Symbol {
    TensorSlotSymmetry::new([(FactorProjector::Antisymmetric, vec![0, 2])])
        .unwrap()
        .register_symbol(
            NamespacedSymbol::parse(name),
            vec![SymbolAttribute::Real],
            vec![],
            UserData::None,
        )
        .unwrap()
}
fn tensor(head: Symbol, a: usize, spectator: usize, b: usize) -> Tensor {
    Tensor::infer(function!(head, slot(a), slot(spectator), slot(b))).unwrap()
}
fn port(tensor: &Tensor, n: usize) -> usize {
    tensor
        .structure()
        .logical_slots()
        .iter()
        .position(|slot| slot.aind == PartialIndex::Explicit(AbstractIndex::Normal(n)))
        .unwrap()
}
fn vector(n: usize) -> Tensor {
    Tensor::infer(function!(
        tensor_symbol!("slot_symmetry_native::v";Real),
        slot(n)
    ))
    .unwrap()
}

#[test]
fn selected_slot_reindexing_preserves_derivative_spectator_and_typed_trace_zero() {
    let head = head("slot_symmetry_native::reindex");
    let value = tensor(head, 1, 7, 2);
    let reversed = value
        .reindex_interface_ports(&HashMap::from([
            (port(&value, 1), AbstractIndex::Normal(2)),
            (port(&value, 2), AbstractIndex::Normal(1)),
        ]))
        .unwrap();
    assert_eq!(value.expression(), &(-reversed.expression()));
    assert_eq!(reversed.rank(), 3);
    let zero = value
        .reindex_interface_ports(&HashMap::from([(
            port(&value, 2),
            AbstractIndex::Normal(1),
        )]))
        .unwrap();
    assert!(zero.is_zero());
    assert_eq!(zero.rank(), 3);
    let finished = zero.clone().finish_parts().unwrap();
    assert!(finished.is_zero());
    assert_eq!(finished.structure(), zero.structure());
    let positions: Vec<_> = zero
        .structure()
        .logical_slots()
        .iter()
        .enumerate()
        .filter_map(|(position, slot)| {
            (slot.aind == PartialIndex::Explicit(AbstractIndex::Normal(1))).then_some(position)
        })
        .collect();
    let trace = zero
        .trace_ports(MatrixChannel {
            input: positions[0],
            output: positions[1],
        })
        .unwrap();
    assert!(trace.is_zero());
    assert_eq!(trace.rank(), 1);
    assert_eq!(
        trace.structure().logical_slots()[0].aind,
        PartialIndex::Explicit(AbstractIndex::Normal(7))
    );
}

#[test]
fn selected_slot_open_occurrences_are_distinct_until_materialization() {
    let head = head("slot_symmetry_native::open");
    let rep = LibraryRep::from(Minkowski {}).new_rep(4).to_symbolic([]);
    let raw = function!(head, rep.clone(), rep.clone(), rep);
    assert!(!raw.is_zero());
    let value = Tensor::infer(raw).unwrap();
    assert_eq!(value.rank(), 3);
    assert!(!value.is_zero());
}

#[test]
fn selected_slot_automorphism_zero_does_not_include_same_representation_spectator() {
    let value = tensor(head("slot_symmetry_native::signed_zero"), 1, 7, 2);
    let zero = value
        .multiply(&vector(1))
        .unwrap()
        .multiply(&vector(2))
        .unwrap()
        .canonize()
        .unwrap();
    assert!(zero.is_zero());
    assert_eq!(zero.rank(), 1);
    let nonzero = value
        .multiply(&vector(1))
        .unwrap()
        .multiply(&vector(7))
        .unwrap()
        .canonize()
        .unwrap();
    assert!(!nonzero.is_zero());
    assert_eq!(nonzero.rank(), 1);
}

#[test]
fn selected_slot_canonicalization_reuses_subgroups_for_dummy_names() {
    let head = head("slot_symmetry_native::dummy");
    let w = |n| {
        Tensor::infer(function!(
            tensor_symbol!("slot_symmetry_native::w";Real),
            slot(n)
        ))
        .unwrap()
    };
    let left = tensor(head, 1, 7, 2)
        .multiply(&vector(1))
        .unwrap()
        .multiply(&w(2))
        .unwrap()
        .canonize()
        .unwrap();
    let right = tensor(head, 8, 7, 9)
        .multiply(&vector(9))
        .unwrap()
        .multiply(&w(8))
        .unwrap()
        .canonize()
        .unwrap();
    assert_eq!(
        (left.expression() + right.expression()).expand(),
        Atom::Zero
    );
    assert_eq!(left.rank(), 1);
    assert_eq!(left.expression(), left.canonize().unwrap().expression());
}

#[test]
fn selected_slot_odd_cycle_handles_repeated_heads_with_different_storage_orders() {
    let head = head("slot_symmetry_native::odd_cycle");
    let value = Tensor::infer(
        function!(head, slot(1), slot(7), slot(2))
            * function!(head, slot(2), slot(3), slot(9))
            * function!(head, slot(9), slot(4), slot(1))
            * vector(7).expression()
            * vector(3).expression()
            * vector(4).expression(),
    )
    .unwrap();
    assert_eq!(value.rank(), 0);
    assert!(value.canonize().unwrap().is_zero());
}

#[test]
fn selected_slot_conjugation_composition_and_atom_reentry_keep_physical_head() {
    let value = tensor(head("slot_symmetry_native::adjoint"), 1, 7, 2);
    assert_eq!(
        value.expression(),
        value.dirac_adjoint(true).unwrap().expression()
    );
    let matrix = Tensor::infer(function!(
        tensor_symbol!("slot_symmetry_native::matrix";Real),
        slot(2),
        slot(3)
    ))
    .unwrap();
    let composed = value
        .compose(
            &matrix,
            MatrixChannel {
                input: port(&value, 1),
                output: port(&value, 2),
            },
            MatrixChannel {
                input: port(&matrix, 2),
                output: port(&matrix, 3),
            },
        )
        .unwrap()
        .canonize()
        .unwrap();
    let product = value.multiply(&matrix).unwrap().canonize().unwrap();
    assert_eq!(composed.expression(), product.expression());
    let admitted = Tensor::infer(composed.expression().clone()).unwrap();
    assert_eq!(admitted.rank(), 3);
    assert_eq!(admitted.expression(), composed.expression());
}

#[test]
fn selected_slot_internal_spectator_dualizes_without_joining_lorentz_group() {
    use crate::representations::ColorFundamental;
    use spenso::structure::representation::Representation;
    let internal: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
    let head = head("slot_symmetry_native::color_spectator");
    let value = Tensor::infer(function!(
        head,
        slot(1),
        internal
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(9))
            .to_atom(),
        slot(2)
    ))
    .unwrap();
    let conjugate = value.dirac_adjoint(true).unwrap();
    assert_eq!(conjugate.rank(), 3);
    assert_eq!(
        conjugate.structure().logical_slots()[port(&conjugate, 9)].rep(),
        internal.dual()
    );
    assert_eq!(
        conjugate.dirac_adjoint(true).unwrap().expression(),
        value.expression()
    );
    let reversed = conjugate
        .reindex_interface_ports(&HashMap::from([
            (port(&conjugate, 1), AbstractIndex::Normal(2)),
            (port(&conjugate, 2), AbstractIndex::Normal(1)),
        ]))
        .unwrap();
    assert_eq!(conjugate.expression(), &(-reversed.expression()));
}
