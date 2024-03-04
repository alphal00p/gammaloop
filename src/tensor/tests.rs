use std::collections::HashMap;

use crate::tensor::{
    ufo::mink_four_vector, Contract, DenseTensor, GetTensorData, HasTensorData, MixedTensor,
    Representation, SparseTensor, SymbolicContract, TensorStructure,
};
use ahash::{HashSet, HashSetExt};
use indexmap::{IndexMap, IndexSet};
use itertools::Itertools;
use num::Complex;
use rand::{rngs, Rng, SeedableRng};
use rand_xoshiro::Xoroshiro64Star;
use smartstring::alias::String;
use symbolica::{
    representations::Atom,
    state::{State, Workspace},
};

use super::{
    structure, symbolic::SymbolicTensor, ufo, AbstractIndex, DataTensor, Dimension,
    HistoryStructure, NumTensor, SetTensorData, Slot, TensorNetwork, VecStructure,
};

fn test_tensor<D, S>(structure: S, seed: u64) -> SparseTensor<D, S>
where
    S: TensorStructure,
    rand::distributions::Standard: rand::distributions::Distribution<D>,
{
    let mut rng: Xoroshiro64Star = Xoroshiro64Star::seed_from_u64(seed);

    let mut tensor = SparseTensor::empty(structure);

    let density = rng.gen_range(0..tensor.size());

    for _ in 0..density {
        tensor
            .set_flat(rng.gen_range(0..tensor.size()), rng.gen())
            .unwrap();
    }

    tensor
}

fn test_structure(length: usize, seed: u64) -> Vec<Slot> {
    let mut rng = Xoroshiro64Star::seed_from_u64(seed);
    let mut s = IndexSet::new();

    let rank = length;
    while s.len() < rank {
        let rep = rng.gen_range(0..=1);
        let dim = Dimension(rng.gen_range(1..=7));
        let id = AbstractIndex(rng.gen_range(0..256));
        let rep = match rep {
            0 => Representation::Euclidean(dim),
            _ => Representation::Lorentz(dim),
        };

        s.insert((id, rep).into());
    }

    s.into_iter().collect_vec()
}

#[test]
fn rng_is_deterministic() {
    let valid = IndexMap::from([
        (vec![1, 0, 1], -7),
        (vec![2, 0, 2], -106),
        (vec![0, 1, 0], -52),
        (vec![2, 1, 2], 4),
        (vec![3, 1, 1], -29),
        (vec![0, 0, 2], 88),
        (vec![2, 0, 1], 79),
        (vec![0, 1, 1], -69),
        (vec![1, 1, 1], -78),
        (vec![3, 0, 2], 55),
        (vec![2, 0, 0], 26),
        (vec![3, 1, 2], -64),
        (vec![1, 0, 2], 68),
    ]);
    for _ in 0..10 {
        let a = test_structure(3, 11);

        let a: SparseTensor<i8> = test_tensor(a, 1);

        assert_eq!(a.hashmap(), valid);
    }
}

#[test]
fn indexflatten() {
    let a = test_structure(4, 31);
    let idx = vec![1, 2, 3, 1];
    let flatidx = a.flat_index(&idx).unwrap();
    assert_eq!(idx, a.expanded_index(flatidx).unwrap());
}

#[test]
fn trace() {
    let structura =
        HistoryStructure::from_integers(&[(1, 5), (1, 5)].map(|(a, d)| (a.into(), d.into())), "a");
    let a = test_tensor::<i8, _>(structura, 3);
    let f = a.internal_contract();

    assert!(f.is_scalar());
    assert_eq!(f.data(), vec![79]);
}

#[test]
fn construct_dense_tensor() {
    let a = test_structure(4, 32);
    let data = vec![1.0; a.size()];
    let tensor = super::DenseTensor::from_data(&data, a).unwrap();
    let num_tensor: NumTensor = tensor.clone().into();
    let data_tensor: DataTensor<f64, _> = tensor.clone().into();
    let mixed_tensor: MixedTensor<_> = tensor.clone().into();

    assert_eq!(mixed_tensor.try_as_float().unwrap().data(), data);
    assert_eq!(data_tensor.data(), data);
    assert_eq!(num_tensor.try_as_float().unwrap().data(), data);
}

#[test]
fn construct_sparse_tensor() -> Result<(), String> {
    let structure = test_structure(3, 11);
    println!("{:?}", structure);

    let mut a = SparseTensor::empty(structure);
    a.set(&[1, 2, 1], 1)?;
    a.set(&[0, 2, 2], 2)?;
    a.set(&[1, 2, 2], 3)?;
    a.set(&[1, 0, 2], 4)?;
    Ok(())
}

#[test]
fn dense_tensor_shape() {
    let a = HistoryStructure::from_integers(
        &[(1, 2), (2, 3), (3, 4)].map(|(a, d)| (a.into(), d.into())),
        "a",
    );

    let data = vec![1.0; a.size()];
    let a = super::DenseTensor::from_data(&data, a).unwrap();
    assert_eq!(a.shape(), vec![2, 3, 4]);
}

#[test]
fn contract_densor() {
    let structur_a = HistoryStructure::new(
        &[
            (3.into(), Representation::Euclidean(2.into())),
            (1.into(), Representation::Euclidean(2.into())),
        ],
        "a",
    );
    let structur_b = HistoryStructure::new(
        &[
            (2.into(), Representation::Euclidean(2.into())),
            (3.into(), Representation::Euclidean(2.into())),
        ],
        "b",
    );

    let a = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_a).unwrap();
    let b = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_b).unwrap();
    let f = a.contract(&b).unwrap();
    assert_eq!(f.data, [7.0, 10.0, 15.0, 22.0]);

    let structur_a =
        HistoryStructure::from_integers(&[(1, 2), (3, 2)].map(|(a, d)| (a.into(), d.into())), "a");
    let structur_b =
        HistoryStructure::from_integers(&[(3, 2), (4, 2)].map(|(a, d)| (a.into(), d.into())), "b");

    let im = Complex::new(0.0, 1.0);
    let re = Complex::new(1.0, 0.0);
    let a = DenseTensor::from_data(&[1.0 * im, 2.0 * im, 3.0 * im, 4.0 * im], structur_a).unwrap();
    let b = DenseTensor::from_data(&[1.0 * im, 2.0 * im, 3.0 * im, 4.0 * im], structur_b).unwrap();
    let f = a.contract(&b).unwrap();

    assert_eq!(f.data, [-7.0 * re, -10.0 * re, -15.0 * re, -22.0 * re]);
}

#[test]
fn multi_index_contract() {
    let structur_a = HistoryStructure::new(
        &[
            (3.into(), Representation::Lorentz(2.into())),
            (1.into(), Representation::Lorentz(2.into())),
        ],
        "a",
    );
    let structur_b = HistoryStructure::new(
        &[
            (1.into(), Representation::Lorentz(2.into())),
            (3.into(), Representation::Lorentz(2.into())),
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
    let structur_a =
        HistoryStructure::from_integers(&[(1, 2), (3, 2)].map(|(a, d)| (a.into(), d.into())), "a");
    let a = DenseTensor::from_data(&[1.0, 2.0, 3.0, 4.0], structur_a).unwrap();
    let b = a.to_sparse();
    let c = b.to_dense();
    assert_eq!(a, c);
}

#[test]
fn gamma() {
    let g1 = ufo::gamma::<f64>(0.into(), (0.into(), 1.into()));
    let g2 = ufo::gamma::<f64>(1.into(), (1.into(), 2.into()));
    let g3 = ufo::gamma::<f64>(2.into(), (2.into(), 0.into()));

    let c = g1.contract(&g2).unwrap().contract(&g3).unwrap();
    println!("Sparse: {:?}", c.data());

    let g1d: NumTensor = g1.to_dense().into();
    let g2d: NumTensor = g2.to_dense().into();
    let g3d: NumTensor = g3.into();

    let cdense: NumTensor = g1d.contract(&g2d).unwrap().contract(&g3d).unwrap();

    println!("{:?}", cdense.try_as_complex().unwrap().data());
    let d = ufo::gamma::<f32>(0.into(), (0.into(), 0.into())).internal_contract();

    println!("{:?}", d.data());
}

#[test]
fn matches() {
    let structur_a = HistoryStructure::new(
        &[
            (3.into(), Representation::Lorentz(2.into())),
            (2.into(), Representation::Lorentz(3.into())),
            (2.into(), Representation::Euclidean(2.into())),
            (1.into(), Representation::Lorentz(2.into())),
        ],
        "a",
    );
    let structur_b = HistoryStructure::new(
        &[
            (1.into(), Representation::Lorentz(2.into())),
            (3.into(), Representation::Lorentz(2.into())),
            (2.into(), Representation::Lorentz(2.into())),
            (1.into(), Representation::Euclidean(2.into())),
        ],
        "b",
    );

    let a = structur_a.match_index(&structur_b);

    println!("{a:?}");
}

#[test]
fn mixed_tensor_contraction() {
    let im = Complex::new(1.5, 1.25);
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let structur_a =
        HistoryStructure::from_integers(&[(2, 2), (1, 2)].map(|(a, d)| (a.into(), d.into())), "a");

    let a = SparseTensor::from_data(&data_a, structur_a.clone()).unwrap();

    let structur_b =
        HistoryStructure::from_integers(&[(2, 2), (4, 2)].map(|(a, d)| (a.into(), d.into())), "b");

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
    let a: NumTensor = ufo::gamma(1.into(), (2.into(), 3.into())).into();
    let b: NumTensor = ufo::gamma(2.into(), (3.into(), 4.into())).into();
    let c: NumTensor = ufo::gamma(3.into(), (4.into(), 2.into())).into();
    let p: NumTensor = mink_four_vector(2.into(), &[2., 3., 2., 1.]).into();
    let q: NumTensor = mink_four_vector(3.into(), &[2., 3., 2., 1.]).into();

    let mut n = TensorNetwork::from(vec![a, b, c, p, q]);

    println!("{}", n.dot());

    if n.graph.validate_neighbors() {
        println!("Validated")
    } else {
        println!("Not validated")
    }

    println!("{:?}", n.graph.neighbors.len());

    n.contract();

    // println!("{}", n.dot());
}

#[test]
fn sparsedensedensesparse() {
    let a: NumTensor = ufo::gamma(1.into(), (2.into(), 3.into())).into();
    let p: NumTensor = mink_four_vector(1.into(), &[2., 3., 2., 1.]).into();

    println!("{:?}", a.contract(&p).unwrap());
    println!("{:?}", p.contract(&a).unwrap());
}

#[test]
fn contract_spensor() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];
    let structur_a =
        HistoryStructure::from_integers(&[(2, 2), (1, 2)].map(|(a, d)| (a.into(), d.into())), "a");

    let a = SparseTensor::from_data(&data_a, structur_a).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_b =
        HistoryStructure::from_integers(&[(1, 2), (3, 2)].map(|(a, d)| (a.into(), d.into())), "b");

    let b = SparseTensor::from_data(&data_b, structur_b).unwrap();

    let _f = a.contract(&b).unwrap();

    let _result = IndexMap::from([(vec![0, 1], 2.0), (vec![1, 0], 2.0)]);

    // assert_eq!(f.elements, result)
}

#[test]
fn sparse_addition() {
    let data_a = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_a =
        HistoryStructure::from_integers(&[(2, 2), (1, 2)].map(|(a, d)| (a.into(), d.into())), "a");

    let a = SparseTensor::from_data(&data_a, structur_a).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_b =
        HistoryStructure::from_integers(&[(1, 2), (2, 2)].map(|(a, d)| (a.into(), d.into())), "b");

    let b = SparseTensor::from_data(&data_b, structur_b).unwrap();

    let _f = a + b;

    let _result = IndexMap::from([(vec![0, 1], 3.0), (vec![1, 0], 3.0)]);

    // assert_eq!(f.elements, result)
}

#[test]
fn sparse_sub() {
    let data_a = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_a =
        HistoryStructure::from_integers(&[(2, 2), (1, 2)].map(|(a, d)| (a.into(), d.into())), "a");

    let a = SparseTensor::from_data(&data_a, structur_a).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 3.0)];

    let structur_b =
        HistoryStructure::from_integers(&[(2, 2), (1, 2)].map(|(a, d)| (a.into(), d.into())), "a");

    let b = SparseTensor::from_data(&data_b, structur_b).unwrap();

    let _f = a - b;

    let _result = IndexMap::from([(vec![0, 1], -1.0), (vec![1, 0], 0.0)]);
    // assert_eq!(f.elements, result);
    // println!("{:?}", f);
}

#[test]
fn contract_densor_with_spensor() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let structur_a =
        HistoryStructure::from_integers(&[(2, 2), (1, 2)].map(|(a, d)| (a.into(), d.into())), "a");

    let a = SparseTensor::from_data(&data_a, structur_a).unwrap();

    let data_b = [1.0, 2.0, 3.0, 4.0];
    let structur_b =
        HistoryStructure::from_integers(&[(1, 2), (4, 2)].map(|(a, d)| (a.into(), d.into())), "b");

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
    let structur_b =
        HistoryStructure::from_integers(&[(1, 2), (4, 3)].map(|(a, d)| (a.into(), d.into())), "b");
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
//     );    /// let state = State:

//     // symb.contract_with_dense(&a);
//     // let structurb = TensorStructure::from_integers(&[2, 4], &[2, 3]);
//     // let b = DenseTensor::symbolic_labels("b", structurb, &ws, &mut state);
// }

#[test]
fn empty_densor() {
    let empty_structure = Vec::<Slot>::new();

    let empty: DenseTensor<f64> = DenseTensor::default(empty_structure);

    assert_eq!(*empty.get(&[]).unwrap(), 0.0);
}

#[test]
fn symbolic_contract() {
    let mut state = State::new();
    let ws = Workspace::new();

    let structura = HistoryStructure::from_integers(
        &[(1, 2), (4, 3)].map(|(a, d)| (a.into(), d.into())),
        "T".to_string(),
    );

    let structurb = HistoryStructure::from_integers(
        &[(3, 2), (2, 3)].map(|(a, d)| (a.into(), d.into())),
        "P".to_string(),
    );

    let a = SymbolicTensor::from_named(&structura, &mut state, &ws).unwrap();
    let b = SymbolicTensor::from_named(&structurb, &mut state, &ws).unwrap();
    let f = a.contract_sym(&b, &state, &ws).unwrap();

    // println!("{:?}", f);
    assert_eq!(
        *f.get_atom(),
        Atom::parse("T(euc(2,1),euc(3,4))*P(euc(2,3),euc(3,2))", &mut state, &ws).unwrap()
    );
}
