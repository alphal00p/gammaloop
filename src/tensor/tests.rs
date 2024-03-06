

use crate::{
    tensor::{
        ufo::mink_four_vector, Contract, DenseTensor, GetTensorData, HasTensorData,
        MixedTensor, Representation, SparseTensor, SymbolicContract, TensorStructure,
    },
};
use ahash::{HashSetExt};

use indexmap::{IndexMap, IndexSet};
use itertools::Itertools;
use num::{Complex};

use rand::{distributions::Uniform, Rng, SeedableRng};
use rand_xoshiro::Xoroshiro64Star;
use smartstring::alias::String;
use symbolica::{
    representations::Atom,
    state::{State, Workspace},
};

use super::{
    symbolic::SymbolicTensor, ufo, AbstractIndex, DataTensor, Dimension,
    HistoryStructure, NumTensor, SetTensorData, Slot,
    TensorNetwork, VecStructure,
};

trait Average {
    fn mean(a: Self, b: Self) -> Self;
}

fn test_tensor<D, S>(structure: S, seed: u64, range: Option<(D, D)>) -> SparseTensor<D, S>
where
    S: TensorStructure,
    D: rand::distributions::uniform::SampleUniform,
    Uniform<D>: Copy,

    rand::distributions::Standard: rand::distributions::Distribution<D>,
{
    let mut rng: Xoroshiro64Star = Xoroshiro64Star::seed_from_u64(seed);

    let mut tensor = SparseTensor::empty(structure);

    let density = rng.gen_range(0..tensor.size());

    if let Some((low, high)) = range {
        let multipliable = Uniform::new(low, high);
        for _ in 0..density {
            tensor
                .set_flat(rng.gen_range(0..tensor.size()), rng.sample(multipliable))
                .unwrap();
        }
    } else {
        for _ in 0..density {
            tensor
                .set_flat(rng.gen_range(0..tensor.size()), rng.gen())
                .unwrap();
        }
    }

    tensor
}

fn test_structure(length: usize, seed: u64) -> Vec<Slot> {
    let mut rng = Xoroshiro64Star::seed_from_u64(seed);
    let mut s = IndexSet::new();

    let rank = length;
    while s.len() < rank {
        let rep = rng.gen_range(0..=1);
        let dim = Dimension(rng.gen_range(1..=9));
        let id = AbstractIndex(rng.gen_range(0..256));
        let rep = match rep {
            0 => Representation::Euclidean(dim),
            _ => Representation::Lorentz(dim),
        };

        s.insert((id, rep).into());
    }

    s.into_iter().collect_vec()
}

fn test_structure_with_dims(dims: &[usize], seed: u64) -> Vec<Slot> {
    let mut s = IndexSet::new();
    let mut rng = Xoroshiro64Star::seed_from_u64(seed);

    for d in dims {
        loop {
            let dim: Dimension = (*d).into();
            let rep = rng.gen_range(0..=1);
            let id = AbstractIndex(rng.gen_range(0..256));

            let rep = match rep {
                0 => Representation::Euclidean(dim),
                _ => Representation::Lorentz(dim),
            };

            if s.insert((id, rep).into()) {
                break;
            }
        }
    }

    s.into_iter().collect_vec()
}

#[test]
fn rng_is_deterministic() {
    let valid = IndexMap::from([
        (vec![3, 0, 3], 53),
        (vec![1, 1, 0], 45),
        (vec![1, 1, 2], -99),
        (vec![2, 1, 0], -59),
        (vec![0, 1, 1], -93),
        (vec![2, 0, 1], 105),
        (vec![4, 1, 0], 125),
        (vec![1, 0, 0], -118),
        (vec![0, 0, 3], -26),
        (vec![4, 0, 0], 59),
        (vec![3, 1, 2], 84),
        (vec![3, 1, 0], -13),
        (vec![1, 0, 3], 119),
        (vec![0, 1, 2], 48),
        (vec![1, 0, 2], 17),
        (vec![0, 0, 0], 34),
        (vec![3, 0, 2], 20),
        (vec![4, 0, 2], -3),
        (vec![3, 1, 3], 69),
        (vec![4, 0, 1], 125),
    ]);
    for _ in 0..10 {
        let a = test_structure(3, 11);

        let a: SparseTensor<i8> = test_tensor(a, 1, None);

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
    let a = test_tensor::<i8, _>(structura, 3, None);
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
    a.set(&[1, 0, 1], 1.)?;
    a.set(&[0, 0, 2], 2.)?;
    a.set(&[1, 1, 2], 3.)?;
    a.set(&[1, 0, 2], 4.)?;

    let num_tensor: NumTensor = a.clone().into();
    let data_tensor: DataTensor<f64, _> = a.clone().into();
    let mixed_tensor: MixedTensor<_> = a.clone().into();

    assert_eq!(
        num_tensor.try_as_float().unwrap().hashmap(),
        data_tensor.hashmap()
    );
    assert_eq!(data_tensor.hashmap(), a.hashmap());
    assert_eq!(mixed_tensor.try_as_float().unwrap().hashmap(), a.hashmap());

    Ok(())
}

#[test]
fn tensor_structure_forwarding() {
    let a = test_structure(6, 1);
    let range = Some((-1000, 1000));

    let sparse: SparseTensor<i16> = test_tensor(a.clone(), 1, range);
    let dense: DenseTensor<i16> = test_tensor(a.clone(), 2, range).to_dense();

    assert_eq!(a.strides(), sparse.strides());
    assert_eq!(dense.reps(), a.reps());
}

#[test]
fn scalar_and_dim1_conract() {
    let common = test_structure_with_dims(&[1, 3, 1, 2], 6);
    let structa = [test_structure(2, 32), common.clone()].concat();
    let structb = [test_structure(3, 22), common.clone()].concat();
    let range = Some((-100, 100));

    let mut tensor_1: SparseTensor<i16> = test_tensor(structa, 3, range);
    tensor_1.set_flat(0, 45);
    let mut tensor_2: SparseTensor<i16> = test_tensor(structb, 2, range);
    tensor_2.set_flat(0, 2);
    let f = tensor_1.contract(&tensor_2).unwrap();

    println!("{:?}", f.hashmap());
}
#[test]
fn contract_with_rank_one_in_middle() {
    let s = 12;
    let common = test_structure(3, s);
    let structa: VecStructure = [test_structure(2, s + 1), common.clone()].concat().into();
    let structb: VecStructure = [test_structure(1, s + 2), common].concat().into();

    println!("seed: {s}");

    println!("--");
    println!("{structa}");
    println!("--");
    println!("{structb}");
    let range = Some((-1000, 1000));
    let tensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa, s + 3, range);
    let dense_a: DenseTensor<i32, VecStructure> = tensor_a.to_dense();
    let tensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb, s + 4, range);
    let dense_b: DenseTensor<i32, VecStructure> = tensor_b.to_dense();

    let f = tensor_b.contract(&tensor_a).unwrap().to_dense();
    let g = dense_b.contract(&dense_a).unwrap();

    println!("--");
    println!("{}", f.structure());
    println!("{:?}", f.strides());

    assert_eq!(f.data, g.data);
}

fn test_structure_with_id<T>(ids: T, seed: u64) -> Vec<Slot>
where
    T: Iterator<Item = usize>,
{
    let mut rng = Xoroshiro64Star::seed_from_u64(seed);
    let mut s = Vec::new();

    for id in ids {
        let rep = rng.gen_range(0..=1);
        let dim = Dimension(rng.gen_range(1..=9));
        let id = AbstractIndex(id);
        let rep = match rep {
            0 => Representation::Euclidean(dim),
            _ => Representation::Lorentz(dim),
        };

        s.push((id, rep).into());
    }
    s
}

#[test]
fn skipping() {
    let dim = Representation::Euclidean(4.into());
    let structa: Vec<Slot> = vec![(1.into(), dim).into(), (2.into(), dim).into()];

    let mut tensora = SparseTensor::empty(structa.clone());
    tensora.set(&[0, 0], 1);
    tensora.set(&[3, 3], 2);
    let structb: Vec<Slot> = vec![(1.into(), dim).into()];
    let tensorb = DenseTensor::from_data(&[1, 2, 3, 4], structb).unwrap();
    let f = tensorb.contract(&tensora).unwrap();
    assert_eq!(f.data, vec![1, 0, 0, 8]);
}

#[test]
fn single_contract() {
    let s = 18;
    let range = Some((-1000, 1000));
    let common = test_structure_with_id(0..1, s);
    let mut structa = test_structure_with_id(1..2, s);
    let mut structb = test_structure_with_id(2..3, s);
    let mut rng = Xoroshiro64Star::seed_from_u64(s);

    structa.insert(rng.gen_range(0..structa.len()), common[0]);
    structb.insert(rng.gen_range(0..structb.len()), common[0]);
    structa.sort();
    let structa: VecStructure = structa.into();
    let structb: VecStructure = structb.into();

    let spensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa.clone(), s + 3, range);

    let densor_a: DenseTensor<i32, VecStructure> = spensor_a.to_dense();
    println!("A={:?}", densor_a);

    let spensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb.clone(), s + 4, range);
    let densor_b: DenseTensor<i32, VecStructure> = spensor_b.to_dense();
    println!("B={:?}", densor_b);

    let dense_dense = densor_b.contract(&densor_a).unwrap();
    println!("A*B {:?}", dense_dense);
    let sparse_sparse = spensor_b.contract(&spensor_a).unwrap().to_dense();
    let dense_sparse = densor_b.contract(&spensor_a).unwrap();
    let sparse_dense = spensor_b.contract(&densor_a).unwrap();

    if dense_dense.data() != sparse_sparse.data() {
        println!("S-S not match at seed: {s}");
    }
    if dense_dense.data() != dense_sparse.data() {
        println!("D-S not match at seed: {s}");
    }
    if dense_dense.data() != sparse_dense.data() {
        println!("S-D not match at seed: {s}");
    }
}

#[test]
fn all_single_contractions() {
    let range = Some((-1000, 1000));

    let mut dseq = vec![];
    let mut sseq = vec![];
    let mut sdeq = vec![];

    for s in 0..10000 {
        let common = test_structure_with_id(0..1, s);
        let mut structa = test_structure_with_id(1..2, s);
        let mut structb = test_structure_with_id(2..3, s);
        let mut rng = Xoroshiro64Star::seed_from_u64(s);

        structa.insert(rng.gen_range(0..structa.len()), common[0]);
        structb.insert(rng.gen_range(0..structb.len()), common[0]);
        structa.sort();
        let structa: VecStructure = structa.into();
        let structb: VecStructure = structb.into();

        let spensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa.clone(), s + 3, range);
        let densor_a: DenseTensor<i32, VecStructure> = spensor_a.to_dense();
        let spensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb.clone(), s + 4, range);
        let densor_b: DenseTensor<i32, VecStructure> = spensor_b.to_dense();

        let dense_dense = densor_b.contract(&densor_a).unwrap();
        // println!("{}", dense_dense.structure());
        let sparse_sparse = spensor_b.contract(&spensor_a).unwrap().to_dense();
        let dense_sparse = densor_b.contract(&spensor_a).unwrap();
        let sparse_dense = spensor_b.contract(&densor_a).unwrap();

        if dense_dense.data() != sparse_sparse.data() {
            sseq.push(s);
            println!("S-S not match at seed: {s}");
        }
        if dense_dense.data() != dense_sparse.data() {
            dseq.push(s);
            println!("D-S not match at seed: {s}");
        }
        if dense_dense.data() != sparse_dense.data() {
            sdeq.push(s);
            println!("S-D not match at seed: {s}");
        }
    }

    println!("S-S:{sseq:?}");
    println!("D-S:{dseq:?}");
    println!("S-D:{sdeq:?}");
}

#[test]
fn all_multi_contractions() {
    let _seeds = [48, 50, 118, 225, 234, 310];
    let range = Some((-1000, 1000));

    let mut dseq = vec![];
    let mut sseq = vec![];
    let mut sdeq = vec![];
    for s in 0..1000 {
        let mut rng = Xoroshiro64Star::seed_from_u64(s);
        let ncommon = rng.gen_range(2..5);

        let common = test_structure_with_id(0..ncommon, s);
        let mut structa = test_structure_with_id(ncommon..ncommon + 1, s);
        let mut structb = test_structure_with_id(ncommon + 1..ncommon + 2, s);

        for c in common {
            structa.insert(rng.gen_range(0..structa.len()), c);
            structb.insert(rng.gen_range(0..structb.len()), c);
        }
        structa.sort();
        let structa: VecStructure = structa.into();
        let structb: VecStructure = structb.into();
        // println!("seed: {s}");

        // println!("seed: {s}");

        // println!("--");
        // println!("{structa}");
        // println!("--");
        // println!("{structb}");
        // println!("--A");
        // println!("{structa}");
        // println!("--B");
        // println!("{structb}");
        // println!("--A*B");
        let spensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa.clone(), s + 3, range);
        let densor_a: DenseTensor<i32, VecStructure> = spensor_a.to_dense();
        let spensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb.clone(), s + 4, range);
        let densor_b: DenseTensor<i32, VecStructure> = spensor_b.to_dense();

        let dense_dense = densor_b.contract(&densor_a).unwrap();
        // println!("{}", dense_dense.structure());
        let sparse_sparse = spensor_b.contract(&spensor_a).unwrap().to_dense();
        let dense_sparse = densor_b.contract(&spensor_a).unwrap();
        let sparse_dense = spensor_b.contract(&densor_a).unwrap();

        if dense_dense.data() != sparse_sparse.data() {
            sseq.push(s);
            println!("S-S not match at seed: {s}");
        }
        if dense_dense.data() != dense_sparse.data() {
            dseq.push(s);
            println!("D-S not match at seed: {s}");
        }
        if dense_dense.data() != sparse_dense.data() {
            sdeq.push(s);
            println!("S-D not match at seed: {s}");
        }
    }

    println!("S-S:{sseq:?}");
    println!("D-S:{dseq:?}");
    println!("S-D:{sdeq:?}");
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
