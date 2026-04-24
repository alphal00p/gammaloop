use crate::{
    Spensor, structure::ConvertibleToSpensoName, structure::SpensoName, to_library_spensor,
};
use spenso::{
    network::{library::symbolic::ExplicitKey, parsing::ShadowedStructure},
    structure::{HasName, HasStructure, ScalarStructure, abstract_index::AbstractIndex},
    tensors::data::{DataTensor, SetTensorData, SparseTensor},
};
use symbolica::symbol;

#[test]
fn tensor_set_name_updates_owned_tensor() {
    let mut tensor = Spensor::from(DataTensor::Sparse(SparseTensor::<
        f64,
        ShadowedStructure<AbstractIndex>,
    >::empty(
        ShadowedStructure::scalar_structure(),
        0.0,
    )));

    tensor.set_name(ConvertibleToSpensoName(SpensoName { name: symbol!("T") }));

    assert_eq!(
        tensor.tensor.structure.name(),
        Some(symbol!("T")),
        "set_name should mutate the owned tensor structure rather than a detached copy",
    );
    assert_eq!(tensor.get_name().map(|name| name.name), Some(symbol!("T")));
}

#[test]
fn tensor_to_library_tensor_preserves_name() {
    let mut tensor = Spensor::from(DataTensor::Sparse(SparseTensor::<
        f64,
        ShadowedStructure<AbstractIndex>,
    >::empty(
        ShadowedStructure::scalar_structure(),
        0.0,
    )));
    tensor.set_name(ConvertibleToSpensoName(SpensoName { name: symbol!("L") }));

    let library_tensor = tensor.to_library_tensor();
    let converted_again = to_library_spensor(&tensor);

    assert_eq!(library_tensor.tensor.structure.name(), Some(symbol!("L")));
    assert_eq!(converted_again.tensor.structure.name(), Some(symbol!("L")));
}

#[test]
fn tensor_density_reports_sparse_fill_ratio() {
    let mut sparse = SparseTensor::<f64, ShadowedStructure<AbstractIndex>>::empty(
        ShadowedStructure::scalar_structure(),
        0.0,
    );
    sparse.set_flat(0.into(), 1.0).unwrap();
    let tensor = Spensor::from(DataTensor::Sparse(sparse));

    assert_eq!(tensor.density().unwrap(), 1.0);
}

#[test]
fn tensor_nnz_reports_sparse_entry_count() {
    let mut sparse = SparseTensor::<f64, ShadowedStructure<AbstractIndex>>::empty(
        ShadowedStructure::scalar_structure(),
        0.0,
    );
    sparse.set_flat(0.into(), 1.0).unwrap();
    let tensor = Spensor::from(DataTensor::Sparse(sparse));

    assert_eq!(tensor.nnz().unwrap(), 1);
}

#[test]
fn library_conversion_changes_structure_key_type() {
    let tensor = Spensor::from(DataTensor::Sparse(SparseTensor::<
        f64,
        ShadowedStructure<AbstractIndex>,
    >::empty(
        ShadowedStructure::scalar_structure(),
        0.0,
    )));

    let library_tensor = tensor.to_library_tensor();
    let _: ExplicitKey<AbstractIndex> = library_tensor.tensor.structure.structure().clone();
}
