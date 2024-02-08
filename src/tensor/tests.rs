use crate::tensor::{
    Contract, DenseTensor, HasTensorData, HasTensorStructure, NumTensor,
    Representation::{self, Lorentz},
    SparseTensor,
};
use indexmap::IndexMap;
use num::Complex;
use smartstring::alias::String;
use symbolica::{
    representations::Atom,
    state::{State, Workspace},
};

use rustc_hash::FxHasher;
type FxIndexMap<K, V> = IndexMap<K, V, FxHasher>;

use super::{
    symbolic_tensor::SymbolicTensor, ufo_spin_tensors, NumTensors, TensorNetwork, TensorSkeleton,
};

fn mink_four_vector<T>(index: usize, p: [T; 4]) -> DenseTensor<T>
where
    T: Clone,
{
    DenseTensor::from_data(
        &p,
        TensorSkeleton::from_idxsing(&[(index, Lorentz(4))], "p".into()),
    )
    .unwrap()
}

#[test]
fn indexflatten() {
    let a = TensorSkeleton::from_integers(&[(1, 3), (2, 4), (3, 5)], "a");

    let idx = vec![1, 2, 3];
    let flatidx = a.flat_index(&idx).unwrap();
    println!("{:?}", a.strides());

    println!("{}", flatidx);
    println!("{:?}", a.expanded_index(flatidx).unwrap());
    assert_eq!(idx, a.expanded_index(flatidx).unwrap());
}

#[test]
fn trace() {
    let structura = TensorSkeleton::from_integers(&[(1, 5), (1, 5)], "a");
    let a = SparseTensor::from_data(
        &[
            (vec![0, 0], 1.0),
            (vec![1, 1], 2.0),
            (vec![1, 2], 4.),
            (vec![2, 1], 3.),
            (vec![2, 2], 5.),
            (vec![3, 3], 6.),
            (vec![4, 4], 7.),
        ],
        structura,
    )
    .unwrap();
    let f = a.internal_contract();

    println!("{:?}", f);
}

#[test]
fn construct_dense_tensor() {
    let a = TensorSkeleton::from_integers(&[(1, 2), (2, 3), (3, 4)], "a");

    let data = vec![1.0; a.size()];
    super::DenseTensor::from_data(&data, a).unwrap();
}

#[test]
fn construct_sparse_tensor() -> Result<(), String> {
    let structure = TensorSkeleton::from_integers(&[(1, 2), (2, 3), (3, 4)], "a".into());

    let mut a: SparseTensor<usize> = SparseTensor::empty(structure);
    a.set(&[1, 2, 1], 1)?;
    a.set(&[0, 2, 3], 2)?;
    a.set(&[1, 2, 3], 3)?;
    a.set(&[1, 0, 3], 4)?;
    Ok(())
}

#[test]
fn dense_tensor_shape() {
    let a = TensorSkeleton::from_integers(&[(1, 2), (2, 3), (3, 4)], "a");

    let data = vec![1.0; a.size()];
    let a = super::DenseTensor::from_data(&data, a).unwrap();
    assert_eq!(a.shape(), vec![2, 3, 4]);
}

#[test]
fn contract_densor() {
    let structur_a = TensorSkeleton::from_idxsing(
        &[
            (3, Representation::Euclidean(2)),
            (1, Representation::Euclidean(2)),
        ],
        "a",
    );
    let structur_b = TensorSkeleton::from_idxsing(
        &[
            (2, Representation::Euclidean(2)),
            (3, Representation::Euclidean(2)),
        ],
        "b",
    );

    let a = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_a).unwrap();
    let b = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_b).unwrap();
    let f = a.contract(&b).unwrap();
    assert_eq!(f.data, [7.0, 10.0, 15.0, 22.0]);

    let structur_a = TensorSkeleton::from_integers(&[(1, 2), (3, 2)], "a");
    let structur_b = TensorSkeleton::from_integers(&[(3, 2), (4, 2)], "b");

    let im = Complex::new(0.0, 1.0);
    let re = Complex::new(1.0, 0.0);
    let a = DenseTensor::from_data(&[1.0 * im, 2.0 * im, 3.0 * im, 4.0 * im], structur_a).unwrap();
    let b = DenseTensor::from_data(&[1.0 * im, 2.0 * im, 3.0 * im, 4.0 * im], structur_b).unwrap();
    let f = a.contract(&b).unwrap();

    assert_eq!(f.data, [-7.0 * re, -10.0 * re, -15.0 * re, -22.0 * re]);
}

#[test]
fn multi_index_contract() {
    let structur_a = TensorSkeleton::from_idxsing(
        &[
            (3, Representation::Lorentz(2)),
            (1, Representation::Lorentz(2)),
        ],
        "a",
    );
    let structur_b = TensorSkeleton::from_idxsing(
        &[
            (1, Representation::Lorentz(2)),
            (3, Representation::Lorentz(2)),
        ],
        "b",
    );

    let a = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_a).unwrap();
    let b = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_b).unwrap();
    let f = a.contract(&b).unwrap();
    println!("{:?}", f.data);
}

#[test]
fn dense_to_sparse() {
    let structur_a = TensorSkeleton::from_integers(&[(1, 2), (3, 2)], "a");
    let a = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_a).unwrap();
    let b = a.to_sparse();
    let c = b.to_dense();
    assert_eq!(a, c);
}

#[test]
fn gamma() {
    let g1 = ufo_spin_tensors::gamma::<f64>(0, (0, 1));
    let g2 = ufo_spin_tensors::gamma::<f64>(1, (1, 2));
    let g3 = ufo_spin_tensors::gamma::<f64>(2, (2, 0));

    let c = g1.contract(&g2).unwrap().contract(&g3).unwrap();
    println!("Sparse: {:?}", c.data());

    let g1d: NumTensors = g1.to_dense().into();
    let g2d: NumTensors = g2.to_dense().into();
    let g3d: NumTensors = g3.into();

    let cdense: NumTensors = g1d.contract(&g2d).unwrap().contract(&g3d).unwrap();

    println!("{:?}", cdense.try_as_complex().unwrap().data());
    let d = ufo_spin_tensors::gamma::<f32>(0, (0, 0)).internal_contract();

    println!("{:?}", d.data());
}

#[test]
fn matches() {
    let structur_a = TensorSkeleton::from_idxsing(
        &[
            (3, Representation::Lorentz(2)),
            (2, Representation::Lorentz(3)),
            (2, Representation::Euclidean(2)),
            (1, Representation::Lorentz(2)),
        ],
        "a",
    );
    let structur_b = TensorSkeleton::from_idxsing(
        &[
            (1, Representation::Lorentz(2)),
            (3, Representation::Lorentz(2)),
            (2, Representation::Lorentz(2)),
            (1, Representation::Euclidean(2)),
        ],
        "b",
    );

    let a = structur_a.match_index(&structur_b);

    println!("{:?}", a);
}

#[test]
fn mixed_tensor_contraction() {
    let im = Complex::new(1.5, 1.25);
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let structur_a = TensorSkeleton::from_integers(&[(2, 2), (1, 2)], "a");

    let a = SparseTensor::from_data(&data_a, structur_a.clone()).unwrap();

    let structur_b = TensorSkeleton::from_integers(&[(2, 2), (4, 2)], "b");

    let b = DenseTensor::from_data(
        &[1.0 * im, 2.0 * im, 3.0 * im, 4.0 * im],
        structur_b.clone(),
    )
    .unwrap();

    let f = b.contract(&a).unwrap();

    assert_eq!(f.data, [1.0 * im, 6.0 * im, 2.0 * im, 8.0 * im]);

    let data_a = [(vec![0, 0], 1.0 * im), (vec![1, 1], 2.0 * im)];

    let a = SparseTensor::from_data(&data_a, structur_a).unwrap();

    let b = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_b).unwrap();

    let f = a.contract(&b).unwrap();
    assert_eq!(f.data, [1.0 * im, 2.0 * im, 6.0 * im, 8.0 * im]);
}

#[test]
fn tensor_net() {
    let a: NumTensors = ufo_spin_tensors::gamma(1, (2, 3)).into();
    let b: NumTensors = ufo_spin_tensors::gamma(2, (3, 4)).into();
    let c: NumTensors = ufo_spin_tensors::gamma(3, (4, 2)).into();
    let p: NumTensors = mink_four_vector(2, [2., 3., 2., 1.]).into();
    let q: NumTensors = mink_four_vector(3, [2., 3., 2., 1.]).into();

    let mut n = TensorNetwork::new(vec![a, b, c, p, q]);

    println!("{}", n.dot());

    n.contract();

    println!("{}", n.dot());
}

#[test]
fn sparsedensedensesparse() {
    let a: NumTensors = ufo_spin_tensors::gamma(1, (2, 3)).into();
    let p: NumTensors = mink_four_vector(1, [2., 3., 2., 1.]).into();

    println!("{:?}", a.contract(&p).unwrap());
    println!("{:?}", p.contract(&a).unwrap());
}

#[test]
fn contract_spensor() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];
    let structur_a = TensorSkeleton::from_integers(&[(2, 2), (1, 2)], "a");

    let a = SparseTensor::from_data(&data_a, structur_a).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_b = TensorSkeleton::from_integers(&[(1, 2), (3, 2)], "b");

    let b = SparseTensor::from_data(&data_b, structur_b).unwrap();

    let f = a.contract(&b).unwrap();

    let result = IndexMap::from([(vec![0, 1], 2.0), (vec![1, 0], 2.0)]);

    // assert_eq!(f.elements, result)
}

#[test]
fn sparse_addition() {
    let data_a = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_a = TensorSkeleton::from_integers(&[(2, 2), (1, 2)], "a");

    let a = SparseTensor::from_data(&data_a, structur_a).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_b = TensorSkeleton::from_integers(&[(1, 2), (2, 2)], "b");

    let b = SparseTensor::from_data(&data_b, structur_b).unwrap();

    let f = a + b;

    let result = IndexMap::from([(vec![0, 1], 3.0), (vec![1, 0], 3.0)]);

    // assert_eq!(f.elements, result)
}

#[test]
fn sparse_sub() {
    let data_a = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_a = TensorSkeleton::from_integers(&[(2, 2), (1, 2)], "a");

    let a = SparseTensor::from_data(&data_a, structur_a).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 3.0)];

    let structur_b = TensorSkeleton::from_integers(&[(2, 2), (1, 2)], "a");

    let b = SparseTensor::from_data(&data_b, structur_b).unwrap();

    let f = a - b;

    let result = IndexMap::from([(vec![0, 1], -1.0), (vec![1, 0], 0.0)]);
    // assert_eq!(f.elements, result);
    // println!("{:?}", f);
}

#[test]
fn contract_densor_with_spensor() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let structur_a = TensorSkeleton::from_integers(&[(2, 2), (1, 2)], "a");

    let a = SparseTensor::from_data(&data_a, structur_a).unwrap();

    let data_b = [1.0, 2.0, 3.0, 4.0];
    let structur_b = TensorSkeleton::from_integers(&[(1, 2), (4, 2)], "b");

    let b = DenseTensor::from_data(&data_b, structur_b).unwrap();

    let f = a.contract(&b).unwrap();

    assert_eq!(f.data, [1.0, 2.0, 6.0, 8.0]);
}

// #[test]
// fn symbolic_zeros() {
//     let mut state = State::new();
//     let ws = Workspace::new();
//     let structure = TensorSkeleton::from_integers(&[(1, 2), (3, 2)], "a");

//     let sym_zeros = DenseTensor::symbolic_zeros(structure.clone());

//     let zeros: DenseTensor<f64> = DenseTensor::default(structure);

//     assert_eq!(sym_zeros, zeros.to_symbolic(&ws, &mut state));
// }

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
    let structur_b = TensorSkeleton::from_integers(&[(1, 2), (4, 3)], "b");
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

// #[test]
// fn symbolic_matrix_mult() {
//     let mut state = State::new();
//     let ws = Workspace::new();

//     let structura = TensorStructure::from_integers(&[1, 4], &[2, 3]);
//     let aatom = DenseTensor::symbolic_labels("a", structura, &ws, &mut state);
//     let structurb = TensorStructure::from_integers(&[4, 1], &[3, 2]);
//     let _batom = DenseTensor::symbolic_labels("b", structurb.clone(), &ws, &mut state);

//     let data_b = [1.5, 2.25, 3.5, -17.125, 5.0, 6.0];
//     let b = DenseTensor::from_data(&data_b, structurb).unwrap();

//     let symb = b.to_symbolic(&ws, &mut state);

//     let f = aatom
//         .builder(&state, &ws)
//         .contract(&symb.builder(&state, &ws));

//     assert_eq!(
//         *f.unwrap().finish().get(&[]).unwrap(),
//         Atom::parse(
//             "3/2*a_0_0+7/2*a_0_1+5*a_0_2+9/4*a_1_0-137/8*a_1_1+6*a_1_2",
//             &mut state,
//             &ws
//         )
//         .unwrap()
//     );

//     // symb.contract_with_dense(&a);
//     // let structurb = TensorStructure::from_integers(&[2, 4], &[2, 3]);
//     // let b = DenseTensor::symbolic_labels("b", structurb, &ws, &mut state);
// }

#[test]
fn empty_densor() {
    let empty_structure = TensorSkeleton::from_integers(&[], "".into());

    let empty: DenseTensor<f64> = DenseTensor::default(empty_structure);

    assert_eq!(*empty.get(&[]).unwrap(), 0.0);
}

#[test]
fn symbolic_contract() {
    let mut state = State::new();
    let ws = Workspace::new();

    let structura = TensorSkeleton::from_integers(&[(1, 2), (4, 3)], "a".into());

    let structurb = TensorSkeleton::from_integers(&[(3, 2), (2, 3)], "b".into());

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
