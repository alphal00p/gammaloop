use std::collections::BTreeMap;

use crate::tensor::{
    ContractableWithDense, ContractableWithSparse, DenseTensor, HasTensorStructure, SparseTensor,
    TensorStructure, VecSlotExtension,
};
use num::Complex;
use num::Rational;
use symbolica::{
    representations::{AsAtomView, Atom},
    state::{ResettableBuffer, State, Workspace},
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

    let structur_a = TensorStructure::from_integers(&[1, 3], &[2, 2]);
    let structur_b = TensorStructure::from_integers(&[3, 4], &[2, 2]);

    let im = Complex::new(0.0, 1.0);
    let re = Complex::new(1.0, 0.0);
    let a = DenseTensor::from_data(&[1.0 * im, 2.0 * im, 3.0 * im, 4.0 * im], structur_a).unwrap();
    let b = DenseTensor::from_data(&[1.0 * im, 2.0 * im, 3.0 * im, 4.0 * im], structur_b).unwrap();
    let f = a.contract_with_dense(&b).unwrap();

    assert_eq!(f.data, [-7.0 * re, -10.0 * re, -15.0 * re, -22.0 * re]);
}
#[test]
fn mixed_tensor_contraction() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let a = SparseTensor::from_data(&data_a, &[2, 1]).unwrap();

    let structur_b = TensorStructure::from_integers(&[2, 4], &[2, 2]);

    let im = Complex::new(1.9, 1.2);

    let b = DenseTensor::from_data(&[1.0 * im, 2.0 * im, 3.0 * im, 4.0 * im], structur_b).unwrap();

    let f = b.contract_with_sparse(&a).unwrap();

    println!("{:?}", f);
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
fn sparse_addition() {
    let data_a = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];

    let a = SparseTensor::from_data(&data_a, &[2, 1]).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];

    let b = SparseTensor::from_data(&data_b, &[1, 2]).unwrap();

    let f = a + b;

    let result = BTreeMap::from([(vec![0, 1], 3.0), (vec![1, 0], 3.0)]);

    assert_eq!(f.elements, result)
}

#[test]
fn sparse_sub() {
    let data_a = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];

    let a = SparseTensor::from_data(&data_a, &[2, 1]).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 3.0)];

    let b = SparseTensor::from_data(&data_b, &[2, 1]).unwrap();

    let f = a - b;
    println!("{:?}", f);
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

#[test]
fn atom_builder() {
    let im = Complex::new(0, 1);
    let mut state = State::new();
    let ws: Workspace = Workspace::new();

    let x = Atom::parse("0", &mut state, &ws).unwrap();
    let y = Atom::parse("y", &mut state, &ws).unwrap();

    let mut xb: symbolica::representations::AtomBuilder<
        '_,
        symbolica::state::BufferHandle<'_, Atom>,
    > = x.builder(&state, &ws);

    xb = (-(xb + &y + &x) * &y * &ws.new_num(3) / &ws.new_num(4)).pow(&ws.new_num(5)) / &y;

    let zero = ws.new_num(0);

    let neutral_summand = zero.builder(&state, &ws);

    // vec![neutral_summand; 5];
    let zero = ws.new_num(0);

    let neutral_summand = zero.builder(&state, &ws);
    let mut result_data = (0..4)
        .map(|_| zero.builder(&state, &ws))
        .collect::<Vec<_>>();

    for a in result_data {
        println!("{}", a.as_atom_view().printer(&state));
    }

    println!("{}", xb.as_atom_view().printer(&state));
}
