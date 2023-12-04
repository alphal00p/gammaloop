use std::collections::BTreeMap;

use crate::tensor::{
    ContractableWithDense, ContractableWithSparse, DenseTensor, HasTensorStructure, SparseTensor,
    TensorStructure, VecSlotExtension,
};

#[test]
fn indexflatten() {
    let a = TensorStructure::from_integers(&[1, 2, 3], &[3, 4, 5]);

    let idx = vec![1, 2, 3];
    let flatidx = a.flat_index(&idx).unwrap();
    println!("{:?}", a.strides());

    println!("{}", flatidx);
    println!("{:?}", a.expanded_index(flatidx).unwrap());
}

#[test]
fn construct_dense_tensor() {
    let a = TensorStructure::from_integers(&[1, 2, 3], &[2, 3, 4]);

    let data = vec![1.0; a.size()];
    let a = super::DenseTensor::from_data(&data, a).unwrap();
    println!("{:?}", a);
}

#[test]
fn construct_sparse_tensor() -> Result<(), String> {
    let structure = TensorStructure::from_integers(&[1, 2, 3], &[2, 3, 4]);

    let mut a: SparseTensor<usize> = SparseTensor::empty(structure);
    a.set(&[1, 2, 1], 1)?;
    a.set(&[0, 2, 3], 2)?;
    a.set(&[1, 2, 3], 3)?;
    a.set(&[1, 0, 3], 4)?;
    println!("{:?}", a);
    Ok(())
}

#[test]
fn dense_tensor_shape() {
    let a = TensorStructure::from_integers(&[1, 2, 3], &[2, 3, 4]);

    let data = vec![1.0; a.size()];
    let a = super::DenseTensor::from_data(&data, a).unwrap();
    assert_eq!(a.shape(), vec![2, 3, 4]);
}

#[test]
fn contract_densor() {
    let structur_a = TensorStructure::from_integers(&[1, 3], &[2, 2]);
    let structur_b = TensorStructure::from_integers(&[3, 4], &[2, 2]);

    let a = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_a).unwrap();
    let b = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_b).unwrap();
    let f = a.contract_with_dense(&b).unwrap();
    assert_eq!(f.data, [7.0, 10.0, 15.0, 22.0]);
}

#[test]
fn contract_spensor() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let a = SparseTensor::from_data(&data_a, &[2, 1]).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];

    let b = SparseTensor::from_data(&data_b, &[1, 3]).unwrap();

    let f = a.contract_with_sparse(&b).unwrap();

    let result = BTreeMap::from([(vec![0, 1], 2.0), (vec![1, 0], 2.0)]);

    assert_eq!(f.elements, result)
}

#[test]
fn contract_densor_with_spensor() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let a = SparseTensor::from_data(&data_a, &[2, 1]).unwrap();
    println!("{:?}", a);

    let data_b = [1.0, 2.0, 3.0, 4.0];
    let structur_b = TensorStructure::from_integers(&[1, 4], &[2, 2]);

    let b = DenseTensor::from_data(&data_b, structur_b).unwrap();
    println!("{:?}", b);
    let f = a.contract_with_dense(&b).unwrap();
    println!("{:?}", f);
}
