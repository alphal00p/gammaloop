use std::collections::BTreeMap;

use crate::tensor::{
    ContractableWithDense, ContractableWithSparse, DenseTensor, HasTensorStructure, SparseTensor,
    TensorStructure, VecSlotExtension,
};
use num::Complex;
use symbolica::{
    representations::{AsAtomView, Atom},
    state::{State, Workspace},
};

use super::symbolic_tensor::SymbolicTensor;

#[test]
fn indexflatten() {
    let a = TensorStructure::from_integers(&[1, 2, 3], &[3, 4, 5]);

    let idx = vec![1, 2, 3];
    let flatidx = a.flat_index(&idx).unwrap();
    // println!("{:?}", a.strides());

    // println!("{}", flatidx);
    // println!("{:?}", a.expanded_index(flatidx).unwrap());
    assert_eq!(idx, a.expanded_index(flatidx).unwrap());
}

#[test]
fn construct_dense_tensor() {
    let a = TensorStructure::from_integers(&[1, 2, 3], &[2, 3, 4]);

    let data = vec![1.0; a.size()];
    let a = super::DenseTensor::from_data(&data, a).unwrap();
}

#[test]
fn construct_sparse_tensor() -> Result<(), String> {
    let structure = TensorStructure::from_integers(&[1, 2, 3], &[2, 3, 4]);

    let mut a: SparseTensor<usize> = SparseTensor::empty(structure);
    a.set(&[1, 2, 1], 1)?;
    a.set(&[0, 2, 3], 2)?;
    a.set(&[1, 2, 3], 3)?;
    a.set(&[1, 0, 3], 4)?;
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

    let im = Complex::new(1.5, 1.25);

    let b = DenseTensor::from_data(&[1.0 * im, 2.0 * im, 3.0 * im, 4.0 * im], structur_b).unwrap();

    let f = b.contract_with_sparse(&a).unwrap();

    assert_eq!(f.data, [1.0 * im, 2.0 * im, 6.0 * im, 8.0 * im]);

    // println!("{:?}", f);
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

    let result = BTreeMap::from([(vec![0, 1], -1.0), (vec![1, 0], 0.0)]);
    assert_eq!(f.elements, result);
    // println!("{:?}", f);
}

#[test]
fn contract_densor_with_spensor() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let a = SparseTensor::from_data(&data_a, &[2, 1]).unwrap();

    let data_b = [1.0, 2.0, 3.0, 4.0];
    let structur_b = TensorStructure::from_integers(&[1, 4], &[2, 2]);

    let b = DenseTensor::from_data(&data_b, structur_b).unwrap();

    let f = a.contract_with_dense(&b).unwrap();

    assert_eq!(f.data, [1.0, 2.0, 6.0, 8.0]);
}

#[test]
fn symbolic_zeros() {
    let mut state = State::new();
    let ws = Workspace::new();
    let structure = TensorStructure::from_integers(&[1, 3], &[2, 2]);

    let sym_zeros = DenseTensor::symbolic_zeros(structure.clone());

    let zeros: DenseTensor<f64> = DenseTensor::default(structure);

    assert_eq!(sym_zeros, zeros.to_symbolic(&ws, &mut state));
}

#[test]
fn convert_sym() {
    let mut state = State::new();
    let ws = Workspace::new();
    let i = Complex::new(0.0, 1.0);
    let mut data_b = vec![i * Complex::from(5.0), Complex::from(2.6) + i];
    data_b.append(
        &mut [3.34, -17.125, 5.0, 6.0]
            .iter()
            .map(|x| Complex::from(*x))
            .collect::<Vec<_>>(),
    );
    let structur_b = TensorStructure::from_integers(&[1, 4], &[2, 3]);
    let b = DenseTensor::from_data(&data_b, structur_b).unwrap();

    let symb = b.to_symbolic(&ws, &mut state);

    let expected_data: Vec<Atom> = [
        "5*𝑖",
        "𝑖+5854679515581645/2251799813685248",
        "940126422213591/281474976710656",
        "-137/8",
        "5",
        "6",
    ]
    .iter()
    .map(|x| Atom::parse(x, &mut state, &ws).unwrap())
    .collect();

    assert_eq!(
        symb.iter().map(|(_, x)| x.clone()).collect::<Vec<_>>(),
        expected_data
    );
}

#[test]
fn symbolic_matrix_mult() {
    let mut state = State::new();
    let ws = Workspace::new();

    let structura = TensorStructure::from_integers(&[1, 4], &[2, 3]);
    let aatom = DenseTensor::symbolic_labels("a", structura, &ws, &mut state);
    let structurb = TensorStructure::from_integers(&[4, 1], &[3, 2]);
    let _batom = DenseTensor::symbolic_labels("b", structurb.clone(), &ws, &mut state);

    let data_b = [1.5, 2.25, 3.5, -17.125, 5.0, 6.0];
    let b = DenseTensor::from_data(&data_b, structurb).unwrap();

    let symb = b.to_symbolic(&ws, &mut state);

    let f = aatom
        .builder(&state, &ws)
        .contract_with_dense(&symb.builder(&state, &ws));

    assert_eq!(
        *f.unwrap().finish().get(&[]).unwrap(),
        Atom::parse(
            "3/2*a_0_0+7/2*a_0_1+5*a_0_2+9/4*a_1_0-137/8*a_1_1+6*a_1_2",
            &mut state,
            &ws
        )
        .unwrap()
    );

    // symb.contract_with_dense(&a);
    // let structurb = TensorStructure::from_integers(&[2, 4], &[2, 3]);
    // let b = DenseTensor::symbolic_labels("b", structurb, &ws, &mut state);
}

#[test]
fn empty_densor() {
    let empty_structure = TensorStructure::from_integers(&[], &[]);

    let empty: DenseTensor<f64> = DenseTensor::default(empty_structure);

    assert_eq!(*empty.get(&[]).unwrap(), 0.0);
}

#[test]
fn symbolic_contract() {
    let mut state = State::new();
    let ws = Workspace::new();

    let structura = TensorStructure::from_integers(&[1, 4], &[2, 3]);

    let structurb = TensorStructure::from_integers(&[3, 2], &[2, 3]);

    let labela = state.get_or_insert_fn("T", None).unwrap();
    let labelb = state.get_or_insert_fn("P", None).unwrap();

    let a = SymbolicTensor::new(structura, labela, &ws, &mut state);
    let b = SymbolicTensor::new(structurb, labelb, &ws, &mut state);
    let f = a.builder(&state, &ws).contract(&b);

    // println!("{:?}", f);
    assert_eq!(
        *f.finish().get_atom(),
        Atom::parse("T(euc(2,1),euc(3,4))*P(euc(2,3),euc(3,2))", &mut state, &ws).unwrap()
    );
}
