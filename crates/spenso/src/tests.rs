#[cfg(feature = "shadowing")]
use ahash::{HashMap, HashMapExt};
use indexmap::{IndexMap, IndexSet};

use log::info;

use insta::{assert_ron_snapshot, assert_snapshot, assert_yaml_snapshot};
use rand::{distributions::Uniform, Rng, SeedableRng};
use rand_xoshiro::Xoroshiro64Star;

#[cfg(feature = "shadowing")]
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    parse, symbol,
};


#[cfg(feature = "shadowing")]
use std::hash::{DefaultHasher, Hash};

use crate::{
    complex::Complex,
    contraction::{Contract, Trace},
    data::{
        DataTensor, DenseTensor, GetTensorData, HasTensorData, NumTensor, SetTensorData,
        SparseTensor,
    },
    iterators::{
        AbstractFiber, CoreExpandedFiberIterator, CoreFlatFiberIterator, Fiber, FiberClass,
        IteratesAlongFibers,
    },
    structure::{
        abstract_index::AbstractIndex,
        concrete_index::{ExpandedIndex, FlatIndex},
        dimension::Dimension,
        representation::{
            Euclidean, LibraryRep, LibrarySlot, Lorentz, Minkowski, RepName, Representation,
        },
        slot::{DualSlotTo, Slot},
        HasStructure, HistoryStructure, NamedStructure, StructureContract, TensorStructure,
        VecStructure,
    },
    upgrading_arithmetic::{FallibleAdd, FallibleAddAssign, FallibleMul, FallibleSub},
};
#[cfg(feature = "shadowing")]
use crate::{
    iterators::IteratableTensor,
    parametric::{MixedTensor, ParamTensor},
    shadowing::test::EXPLICIT_TENSOR_MAP,
    structure::{TensorShell, ToSymbolic},
    symbolic::SymbolicTensor,
    upgrading_arithmetic::TryIntoUpgrade,
};

fn test_tensor<D, S>(structure: S, seed: u64, range: Option<(D, D)>) -> SparseTensor<D, S>
where
    S: TensorStructure,
    D: rand::distributions::uniform::SampleUniform,
    Uniform<D>: Copy,

    rand::distributions::Standard: rand::distributions::Distribution<D>,
{
    let mut rng: Xoroshiro64Star = Xoroshiro64Star::seed_from_u64(seed);

    let mut tensor = SparseTensor::empty(structure);

    let density = rng.gen_range(0..tensor.size().unwrap());

    if let Some((low, high)) = range {
        let multipliable = Uniform::new(low, high);
        for _ in 0..density {
            tensor
                .set_flat(
                    rng.gen_range(0..tensor.size().unwrap()).into(),
                    rng.sample(multipliable),
                )
                .unwrap();
        }
    } else {
        for _ in 0..density {
            tensor
                .set_flat(rng.gen_range(0..tensor.size().unwrap()).into(), rng.gen())
                .unwrap();
        }
    }

    tensor
}

fn test_structure(length: usize, seed: u64) -> VecStructure {
    let mut rng = Xoroshiro64Star::seed_from_u64(seed);
    let mut s = IndexSet::new();

    let rank = length;
    while s.len() < rank {
        let rep = rng.gen_range(0..=1);
        let dim = Dimension::Concrete(rng.gen_range(1..=9));
        let id = AbstractIndex::from(rng.gen_range(0..256));
        let rep: Representation<LibraryRep> = match rep {
            0 => LibraryRep::from(Euclidean {}).rep(dim).cast(),
            _ => LibraryRep::from(Minkowski {}).rep(dim).cast(),
        };

        s.insert(rep.slot(id));
    }

    s.into_iter().collect()
}

fn test_structure_with_dims(dims: &[isize], seed: u64) -> VecStructure {
    let mut s = IndexSet::new();
    let mut rng = Xoroshiro64Star::seed_from_u64(seed);

    for d in dims {
        loop {
            if *d < 0 {
                let dim: Dimension = (-*d as usize).into();
                let rep = rng.gen_range(0..=1);
                let id = AbstractIndex::from(rng.gen_range(0..256));

                let rep: Representation<LibraryRep> = match rep {
                    0 => LibraryRep::from(Euclidean {}).rep(dim).cast(),
                    _ => LibraryRep::from(Lorentz {}).rep(dim).cast(),
                };
                if s.insert(rep.slot(id)) {
                    break;
                }
            } else {
                let dim: Dimension = (*d as usize).into();
                let rep = rng.gen_range(0..=1);
                let id = AbstractIndex::from(rng.gen_range(0..256));

                let rep: Representation<LibraryRep> = match rep {
                    0 => LibraryRep::from(Euclidean {}).rep(dim).cast(),
                    _ => Lorentz {}.dual().rep(dim).cast(),
                };
                if s.insert(rep.slot(id)) {
                    break;
                }
            }
        }
    }

    s.into_iter().collect()
}

#[test]
fn rng_is_deterministic() {
    let a = test_structure(3, 11);

    let a: SparseTensor<i8> = test_tensor(a, 1, None);

    assert_ron_snapshot!(a.data());
    for _ in 0..10 {
        let b = test_structure(3, 11);

        let b: SparseTensor<i8> = test_tensor(b, 1, None);
        assert_eq!(a.data(), b.data());
    }
}

#[test]
fn indexflatten() {
    let a = test_structure(4, 32);
    // println!("{}", a);
    let idx = vec![1, 5, 0, 2];
    let flatidx = a.flat_index(&idx).unwrap();
    assert_eq!(ExpandedIndex::from(idx), a.expanded_index(flatidx).unwrap());
}

#[test]
fn single_fiber() {
    let a = test_structure(5, 5);

    let fiber = Fiber::from_filter(&[true, false, false, false, false], &a);
    assert!(fiber.single().is_some());

    let iter = fiber.clone().iter();

    let fciter: Vec<FlatIndex> = iter.collect();

    let fiberclass: Vec<FlatIndex> = FiberClass::from(fiber)
        .iter()
        .map(|fc| fc.iter.zero_index)
        .collect::<Vec<_>>();

    assert_ron_snapshot!((a, fciter, fiberclass));
}

#[test]
fn expanded_vs_flat_iter() {
    let a = test_structure(5, 5);
    let fiber = Fiber::from_filter(&[false, true, true, false, true], &a);

    let flat_iter = CoreFlatFiberIterator::new(&fiber, false);
    let expanded_iter = CoreExpandedFiberIterator::new(&fiber, false);

    let flat: Vec<FlatIndex> = flat_iter.collect();
    let expanded: Vec<FlatIndex> = expanded_iter.collect();

    assert_eq!(flat, expanded);
}

#[test]
fn fiber_class_vs_fiber_iterator() {
    let a = test_structure(5, 5);
    let fiber_class: FiberClass<'_, VecStructure> =
        Fiber::from_filter(&[false, true, true, false, true], &a).into();

    let fiber = Fiber::from_filter(&[true, false, false, true, false], &a);

    let flat_iter = CoreFlatFiberIterator::new(&fiber, false);
    let flat_class_iter = CoreFlatFiberIterator::new(&fiber_class, false);
    let flat: Vec<FlatIndex> = flat_iter.collect();
    let flat_class: Vec<FlatIndex> = flat_class_iter.collect();

    assert_eq!(flat, flat_class);
    let expanded_iter = CoreExpandedFiberIterator::new(&fiber, false);
    let expanded_class_iter = CoreExpandedFiberIterator::new(&fiber_class, false);

    let expanded: Vec<FlatIndex> = expanded_iter.collect();
    let expanded_class: Vec<FlatIndex> = expanded_class_iter.collect();
    assert_eq!(expanded, expanded_class);
}

#[test]
fn fibers() {
    let a = test_structure(5, 5);
    let fiber = Fiber::from_filter(&[true, true, false, false, true], &a);
    let fiberclass: FiberClass<'_, VecStructure> = fiber.into();
    let iter = fiberclass.clone().iter();
    let fciter: Vec<FlatIndex> = iter.map(|fc| fc.iter.zero_index).collect();

    let fiter: Vec<FlatIndex> = fiberclass
        .iter()
        .flat_map(|fc| fc.collect::<Vec<FlatIndex>>())
        .collect();
    assert_ron_snapshot!((a, fiter, fciter));
}

#[test]
fn fiber_from_structure() {
    let a = VecStructure::from_iter([
        LibraryRep::slot(Euclidean {}.into(), 4, 104),
        Euclidean {}.slot(1, 128).cast(),
        Euclidean {}.slot(5, 164).cast(),
        Lorentz {}.slot(7, 88).cast(),
        Lorentz {}.dual().slot(2, 145).cast(),
    ]);

    let fiber = Fiber::from(([1u8, 1, 0, 1, 0].as_slice()).into(), &a);

    let fiber_iter = fiber
        .clone()
        .iter()
        .map(|f| a.expanded_index(f).unwrap())
        .collect::<Vec<_>>();

    let fiberclass: FiberClass<'_, VecStructure> = fiber.into();

    let fiberclass_iter: Vec<ExpandedIndex> = fiberclass
        .iter()
        .flat_map(|i| i.map(|f| a.expanded_index(f).unwrap()).collect::<Vec<_>>())
        .collect();

    assert_eq!(fiberclass_iter.len(), a.size().unwrap());
    assert_yaml_snapshot!((fiberclass_iter, fiber_iter));
}

#[test]
fn permutation() {
    let a: Vec<LibrarySlot> = vec![
        Euclidean {}.slot(1, 2).cast(),
        Lorentz {}.slot(3, 4).cast(),
        Euclidean {}.slot(5, 6).cast(),
    ];

    let b: Vec<LibrarySlot> = vec![
        Lorentz {}.slot(3, 4).cast(),
        Euclidean {}.slot(5, 6).cast(),
        Euclidean {}.slot(1, 2).cast(),
    ];

    let permutation = a.find_permutation(&b).unwrap();
    // println!("{:?}", permutation);

    let c = permutation.apply_slice(&b);

    assert_eq!(c, a);
}

#[test]
fn trace() {
    let structura: HistoryStructure<String, ()> =
        HistoryStructure::from(NamedStructure::from_iter(
            [Euclidean {}.slot(5, 1), Euclidean {}.slot(5, 1)],
            "a".into(),
            None,
        ));
    let a = test_tensor::<i8, _>(structura, 3, None);
    let f = a.internal_contract();

    assert!(f.is_scalar());
    assert_eq!(f.data(), vec![79]);
}

#[test]
fn construct_dense_tensor() {
    let a = test_structure(4, 32);
    let data = vec![1.0; a.size().unwrap()];
    let tensor = DenseTensor::from_data(data.clone(), a).unwrap();
    let num_tensor: NumTensor = tensor.clone().into();
    let data_tensor: DataTensor<f64, _> = tensor.clone().into();

    #[cfg(feature = "shadowing")]
    {
        let mixed_tensor: MixedTensor<f64, _> = tensor.clone().into();
        assert_eq!(
            mixed_tensor
                .try_into_concrete()
                .unwrap()
                .try_as_real()
                .unwrap()
                .data(),
            data
        );
    }

    assert_eq!(data_tensor.data(), data);
    assert_eq!(num_tensor.try_as_float().unwrap().data(), data);
}

use eyre::Result;
#[test]
fn construct_sparse_tensor() -> Result<()> {
    let structure = test_structure(3, 11);
    // println!("{}", structure);

    let mut a = SparseTensor::empty(structure);
    a.set(&[1, 2, 1], 1.)?;
    a.set(&[0, 0, 0], 2.)?;
    a.set(&[1, 1, 1], 3.)?;
    a.set(&[1, 2, 1], 4.)?;

    let num_tensor: NumTensor = a.clone().into();
    let data_tensor: DataTensor<f64, _> = a.clone().into();

    #[cfg(feature = "shadowing")]
    {
        let mixed_tensor: MixedTensor<f64, _> = a.clone().into();
        assert_eq!(
            mixed_tensor
                .try_into_concrete()
                .unwrap()
                .try_as_real()
                .unwrap()
                .hashmap(),
            a.hashmap()
        );
    }
    assert_eq!(
        num_tensor.try_as_float().unwrap().hashmap(),
        data_tensor.hashmap()
    );
    assert_eq!(data_tensor.hashmap(), a.hashmap());

    Ok(())
}

#[test]
fn tensor_structure_forwarding() {
    let a = test_structure(6, 1);
    let range = Some((-1000, 1000));

    let sparse: SparseTensor<i16> = test_tensor(a.clone(), 1, range);
    let dense: DenseTensor<i16> = test_tensor(a.clone(), 2, range).to_dense();

    assert_eq!(a.strides().unwrap(), sparse.strides().unwrap());
    assert_eq!(dense.reps(), a.reps());
}

#[test]
fn scalar_and_dim1_conract() {
    let common = test_structure_with_dims(&[1, 3, 1, 2], 6);
    let commondual = test_structure_with_dims(&[-1, -3, -1, -2], 6);
    assert_snapshot!(common.to_string());
    let mut structa = test_structure(1, 32);
    assert_snapshot!(structa.to_string());
    structa.merge(&common);
    let mut structb = test_structure(1, 22);
    assert_snapshot!(structb.to_string());
    structb.merge(&commondual);
    let range = Some((-100, 100));

    let mut tensor_1: SparseTensor<i16> = test_tensor(structa, 3, range);
    tensor_1.set_flat(0.into(), 45).unwrap();

    let mut tensor_2: SparseTensor<i16> = test_tensor(structb, 2, range);
    tensor_2.set_flat(0.into(), 2).unwrap();
    let f = tensor_1.contract(&tensor_2).unwrap();

    let g = tensor_1.contract(&tensor_2.to_dense()).unwrap();
    let h = tensor_1.to_dense().contract(&tensor_2).unwrap();

    let i = tensor_1.to_dense().contract(&tensor_2.to_dense()).unwrap();
    assert_eq!(f.to_dense().data(), h.data());
    assert_eq!(f.to_dense().data(), g.data());
    assert_eq!(f.to_dense().data(), i.data());
    assert_snapshot!(f.structure().to_string());

    assert_ron_snapshot!(f.data());
}

#[test]
fn dense_dense_single_contract() {
    let structura: VecStructure<Euclidean> =
        VecStructure::new(vec![Euclidean {}.slot(4, 1), Euclidean {}.slot(4, 2)]);
    let structurb = VecStructure::new(vec![Euclidean {}.slot(4, 2), Euclidean {}.slot(3, 3)]);

    let a = DenseTensor::from_data(
        vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
        structura,
    )
    .unwrap();

    let b = DenseTensor::from_data(
        vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        structurb.clone(),
    )
    .unwrap();

    let f = a.contract(&b).unwrap();

    assert_yaml_snapshot!(f);
}

#[test]
fn sparse_diag_dense_contract() {
    // Logger::try_with_str("trace").unwrap().start().unwrap();
    let structura: VecStructure<Euclidean> =
        VecStructure::new(vec![Euclidean {}.slot(4, 1), Euclidean {}.slot(4, 2)]);
    let structurb = VecStructure::new(vec![Euclidean {}.slot(4, 2), Euclidean {}.slot(3, 3)]);

    let a = SparseTensor::from_data(
        [
            (vec![0, 0], 1),
            (vec![1, 1], 2),
            (vec![2, 2], 3),
            (vec![3, 3], 4),
        ],
        structura.clone(),
    )
    .unwrap();

    let b = DenseTensor::from_data(
        vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        structurb.clone(),
    )
    .unwrap();

    let f = a.contract(&b).unwrap();
    let g = b.contract(&a).unwrap();
    let h = a.contract(&b.to_sparse()).unwrap();
    let i = b.contract(&a.to_dense()).unwrap();

    assert_eq!(i.data(), g.data());
    assert_eq!(i.data(), h.to_dense().data());
    assert_eq!(f.data(), i.data());

    assert_yaml_snapshot!(f);
}

#[test]
fn contract_with_rank_one_in_middle() {
    let s = 12;
    let common = test_structure(3, s);
    let mut structa: VecStructure = test_structure(2, s + 1);
    structa.merge(&common);
    let mut structb: VecStructure = test_structure(1, s + 2);
    structb.merge(&common);

    // println!("seed: {s}");

    // println!("--");
    // println!("{structa}");
    // println!("--");
    // println!("{structb}");
    let range = Some((-1000, 1000));
    let tensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa, s + 3, range);
    let dense_a: DenseTensor<i32, VecStructure> = tensor_a.to_dense();
    let tensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb, s + 4, range);
    let dense_b: DenseTensor<i32, VecStructure> = tensor_b.to_dense();

    let f = tensor_b.contract(&tensor_a).unwrap().to_dense();
    let g = dense_b.contract(&dense_a).unwrap();

    assert_eq!(f.data, g.data);
}

fn test_structure_with_id<T>(ids: T, seed: u64) -> Vec<LibrarySlot>
where
    T: IntoIterator<Item = i32>,
{
    let mut rng = Xoroshiro64Star::seed_from_u64(seed);
    let mut s = Vec::new();

    for id in ids {
        let rep = rng.gen_range(0..=1);
        let dim = Dimension::Concrete(rng.gen_range(1..=9));
        let (rep, id) = if id > 0 {
            let id = AbstractIndex::from(id);
            let rep: Representation<LibraryRep> = match rep {
                0 => Euclidean {}.rep(dim).cast(),
                _ => Lorentz {}.rep(dim).cast(),
            };
            (rep, id)
        } else {
            let id = AbstractIndex::from(-id);
            let rep: Representation<LibraryRep> = match rep {
                0 => Euclidean {}.rep(dim).cast(),
                _ => Lorentz {}.dual().rep(dim).cast(),
            };
            (rep, id)
        };

        s.push(rep.slot(id));
    }
    s
}

#[test]
fn single_contract() {
    // Logger::try_with_str("trace").unwrap().start().unwrap();
    let s = 18;
    let range = Some((-1000, 1000));
    let common = test_structure_with_id(1..2, s);
    let commondual = test_structure_with_id([-1], s);
    let mut structa = test_structure_with_id(2..3, s);
    let mut structb = test_structure_with_id([-3, 4], s);
    let mut rng = Xoroshiro64Star::seed_from_u64(s);

    structa.insert(rng.gen_range(0..structa.len()), common[0]);
    structb.insert(rng.gen_range(0..structb.len()), commondual[0]);
    structa.sort();
    let structa: VecStructure = structa.into();

    let structb: VecStructure = structb.into();
    let spensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa.clone(), s + 3, range);

    let densor_a: DenseTensor<i32, VecStructure> = spensor_a.to_dense();
    // println!("A={:?}", densor_a);

    let spensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb.clone(), s + 4, range);
    assert_ron_snapshot!((structa, structb));
    let densor_b: DenseTensor<i32, VecStructure> = spensor_b.to_dense();
    // println!("B={:?}", densor_b);

    let dense_dense = densor_b.contract(&densor_a).unwrap();
    assert_ron_snapshot!(dense_dense.structure());
    // println!("A*B {:?}", dense_dense);
    let sparse_sparse = spensor_b.contract(&spensor_a).unwrap().to_dense();
    let dense_sparse = densor_b.contract(&spensor_a).unwrap();
    let sparse_dense = spensor_b.contract(&densor_a).unwrap();

    assert_eq!(
        dense_dense.data(),
        sparse_sparse.data(),
        "S-S not match at seed: {s}"
    );
    assert_eq!(
        dense_dense.data(),
        dense_sparse.data(),
        "D-S not match at seed: {s}"
    );
    assert_eq!(
        dense_dense.data(),
        sparse_dense.data(),
        "S-D not match at seed: {s}"
    );
}

#[test]
fn all_single_contractions() {
    // Logger::try_with_str("trace").unwrap().start().unwrap();
    let range = Some((-1000, 1000));

    let mut dseq = vec![];
    let mut sseq = vec![];
    let mut sdeq = vec![];

    for s in 0..1000 {
        let common = test_structure_with_id(1..2, s);
        let commondual = test_structure_with_id([-1], s);
        let mut structa = test_structure_with_id(2..3, s);
        let mut structb = test_structure_with_id([-3, 4], s);
        let mut rng = Xoroshiro64Star::seed_from_u64(s);

        structa.insert(rng.gen_range(0..structa.len()), common[0]);
        structb.insert(rng.gen_range(0..structb.len()), commondual[0]);
        structa.sort();
        let structa: VecStructure = structa.into();
        let structb: VecStructure = structb.into();

        let spensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa.clone(), s + 3, range);
        let densor_a: DenseTensor<i32, VecStructure> = spensor_a.to_dense();
        let spensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb.clone(), s + 4, range);
        let densor_b: DenseTensor<i32, VecStructure> = spensor_b.to_dense();

        let dense_dense = densor_b.contract(&densor_a).unwrap();
        // println!("{}", dense_dense.structure());
        let zeros: DenseTensor<i32> = DenseTensor::zero(dense_dense.structure.clone());
        let sparse_sparse = spensor_b.contract(&spensor_a).unwrap().to_dense();
        let dense_sparse = densor_b.contract(&spensor_a).unwrap_or(zeros.clone());
        let sparse_dense = spensor_b.contract(&densor_a).unwrap_or(zeros.clone());

        if dense_dense.data() != sparse_sparse.data() {
            sseq.push(s);
        }
        if dense_dense.data() != dense_sparse.data() {
            dseq.push(s);
        }
        if dense_dense.data() != sparse_dense.data() {
            sdeq.push(s);
        }
    }

    assert_eq!(sseq.len(), 0, "Sparse-Sparse failed at seeds {sseq:?}");
    assert_eq!(dseq.len(), 0, "Dense-Sparse failed at seeds {dseq:?}");
    assert_eq!(sdeq.len(), 0, "Sparse-Dense failed at seeds {sdeq:?}");
}

#[test]
fn simple_multi_contract() {
    // Logger::try_with_str("trace").unwrap().start().unwrap();
    let structa: VecStructure = VecStructure::<LibraryRep>::from_iter([
        LibraryRep::from(Euclidean {}).slot(3, 1),
        Euclidean {}.slot(4, 2).cast(),
        Lorentz {}.slot(4, 3).cast(),
    ]);
    let structb: VecStructure = VecStructure::new(vec![
        Euclidean {}.slot(4, 2).cast(),
        Euclidean {}.slot(3, 4).cast(),
        Lorentz {}.dual().slot(4, 3).cast(),
    ]);

    let a = DenseTensor::from_data(
        vec![
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 12, 13, 14, 15, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11,
            12, 13, 14, 15, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 12, 13, 14, 15, 16,
        ],
        structa,
    )
    .unwrap();

    let b = DenseTensor::from_data(
        vec![
            3, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 3, 2, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0,
            3, 2, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 3, 2, 1,
        ],
        structb,
    )
    .unwrap();

    let f = a.contract(&b).unwrap();

    let g = b.contract(&a).unwrap();
    assert_eq!(f.data(), g.data(), "A-B not equal to B-A");

    let h = a.contract(&b.to_sparse()).unwrap();
    assert_eq!(f.data(), h.data(), "dense dense not equal to dense sparse");

    let i = b.contract(&a.to_sparse()).unwrap();

    assert_eq!(f.data(), i.data(), "dense dense not equal to sparse dense");
    let j = a.to_sparse().contract(&b.to_sparse()).unwrap();

    assert_eq!(
        f.data(),
        j.to_dense().data(),
        "dense dense not equal to sparse sparse"
    );

    assert_yaml_snapshot!(f.data());
}

#[test]
fn multi_contract_permuted() {
    // Logger::try_with_str("trace").unwrap().start().unwrap();
    let range = Some((-1000, 1000));
    let s = 18;
    let mut rng = Xoroshiro64Star::seed_from_u64(s);
    let ncommon = rng.gen_range(2..5);

    let common = test_structure_with_id(1..ncommon + 1, s);

    let dualcommon = test_structure_with_id((-ncommon..=-1).rev(), s);
    let mut structa = test_structure_with_id(ncommon + 1..ncommon + 2, s);
    let mut structb = test_structure_with_id(ncommon + 2..ncommon + 3, s);

    for (c, dualc) in common.into_iter().zip(dualcommon.into_iter()) {
        structa.insert(rng.gen_range(0..structa.len()), c);
        info!("inserting {c} into structa");
        structb.insert(rng.gen_range(0..structb.len()), dualc);
        info!("inserting {dualc} into structb");
    }
    structa.sort();
    let structa: VecStructure = structa.into();
    let structb: VecStructure = structb.into();

    println!("structa = \n{}", structa);
    println!("structb = \n{}", structb);

    let spensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa.clone(), s + 3, range);
    let densor_a: DenseTensor<i32, VecStructure> = spensor_a.to_dense();
    let spensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb.clone(), s + 4, range);
    let densor_b: DenseTensor<i32, VecStructure> = spensor_b.to_dense();

    let dense_dense = densor_b.contract(&densor_a).unwrap();
    // println!("{}", dense_dense.structure());
    let sparse_sparse = spensor_b.contract(&spensor_a).unwrap().to_dense();
    let dense_sparse = densor_b.contract(&spensor_a).unwrap();
    let sparse_dense = spensor_b.contract(&densor_a).unwrap();

    assert_eq!(
        dense_dense.data(),
        sparse_dense.data(),
        "S-D not match at seed: {s}"
    );
    assert_eq!(
        dense_dense.data(),
        dense_sparse.data(),
        "D-S not match at seed: {s}"
    );
    assert_eq!(
        dense_dense.data(),
        sparse_sparse.data(),
        "S-S not match at seed: {s}"
    );
}

#[test]
fn multi_contract() {
    let range = Some((-1000, 1000));
    let s = 18;
    let mut rng = Xoroshiro64Star::seed_from_u64(s);
    let ncommon = rng.gen_range(2..5);
    let common = test_structure_with_id(1..ncommon + 1, s);

    let dualcommon = test_structure_with_id((-ncommon..=-1).rev(), s);
    let mut structa = test_structure_with_id(ncommon + 1..ncommon + 2, s);
    let mut structb = test_structure_with_id(ncommon + 2..ncommon + 3, s);

    for (c, dualc) in common.into_iter().zip(dualcommon.into_iter()) {
        structa.insert(rng.gen_range(0..structa.len()), c);

        structb.insert(rng.gen_range(0..structb.len()), dualc);
    }

    let structa: VecStructure = structa.into();
    let structb: VecStructure = structb.into();

    println!("structa = {}", structa);
    println!("structb = {}", structb);

    let spensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa.clone(), s + 3, range);
    // insta::assert_ron_snapshot!(spensor_a);
    let densor_a: DenseTensor<i32, VecStructure> = spensor_a.to_dense();
    let spensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb.clone(), s + 4, range);
    // insta::assert_ron_snapshot!(spensor_b);
    let densor_b: DenseTensor<i32, VecStructure> = spensor_b.to_dense();

    let dense_dense = densor_b.contract(&densor_a).unwrap();
    // println!("{}", dense_dense.structure());
    let sparse_sparse = spensor_b.contract(&spensor_a).unwrap().to_dense();
    let dense_sparse = densor_b.contract(&spensor_a).unwrap();
    let sparse_dense = spensor_b.contract(&densor_a).unwrap();

    assert_eq!(
        dense_dense.data(),
        sparse_dense.data(),
        "S-D not match at seed: {s}"
    );
    assert_eq!(
        dense_dense.data(),
        sparse_sparse.data(),
        "S-S not match at seed: {s}"
    );
    assert_eq!(
        dense_dense.data(),
        dense_sparse.data(),
        "D-S not match at seed: {s}"
    );

    insta::assert_ron_snapshot!(dense_dense);
}

#[test]
fn all_multi_contractions() {
    // Logger::try_with_str("trace").unwrap().start().unwrap();
    let _seeds = [48, 50, 118, 225, 234, 310];
    let range = Some((-1000, 1000));

    let mut dseq = vec![];
    let mut sseq = vec![];
    let mut sdeq = vec![];
    for s in 0..1000 {
        let mut rng = Xoroshiro64Star::seed_from_u64(s);
        // let ncommon = rng.gen_range(2..5);
        let ncommon = rng.gen_range(2..5);
        let common = test_structure_with_id(1..ncommon + 1, s);

        let dualcommon = test_structure_with_id((-ncommon..=-1).rev(), s);
        let mut structa = test_structure_with_id(ncommon + 1..ncommon + 2, s);
        let mut structb = test_structure_with_id(ncommon + 2..ncommon + 3, s);

        for (c, dualc) in common.into_iter().zip(dualcommon.into_iter()) {
            structa.insert(rng.gen_range(0..structa.len()), c);

            structb.insert(rng.gen_range(0..structb.len()), dualc);
        }
        structa.sort();
        let structa: VecStructure = structa.into();
        let structb: VecStructure = structb.into();

        let spensor_a: SparseTensor<i32, VecStructure> = test_tensor(structa.clone(), s + 3, range);
        let densor_a: DenseTensor<i32, VecStructure> = spensor_a.to_dense();
        let spensor_b: SparseTensor<i32, VecStructure> = test_tensor(structb.clone(), s + 4, range);
        let densor_b: DenseTensor<i32, VecStructure> = spensor_b.to_dense();

        let dense_dense = densor_b.contract(&densor_a).unwrap();
        // println!("{}", dense_dense.structure());
        let zeros: DenseTensor<i32> = DenseTensor::zero(dense_dense.structure.clone());
        let sparse_sparse = spensor_b.contract(&spensor_a).unwrap().to_dense();
        let dense_sparse = densor_b.contract(&spensor_a).unwrap_or(zeros.clone());
        let sparse_dense = spensor_b.contract(&densor_a).unwrap_or(zeros.clone());

        if dense_dense.data() != sparse_sparse.data() {
            sseq.push(s);
        }
        if dense_dense.data() != dense_sparse.data() {
            dseq.push(s);
        }
        if dense_dense.data() != sparse_dense.data() {
            sdeq.push(s);
        }
    }
    assert_eq!(sdeq.len(), 0, "Sparse-Dense failed at seeds {sdeq:?}");
    // assert_eq!(dseq.len(), 0, "Dense-Sparse failed at seeds {dseq:?}");
    // assert_eq!(sseq.len(), 0, "Sparse-Sparse failed at seeds {sseq:?}");
}

// #[test]
// fn gamma() {
//     let g1: SparseTensor<Complex<f64>> = ufo::gamma(0.into(), [0.into(), 1.into()]);
//     let g2: SparseTensor<Complex<f64>> = ufo::gamma(1.into(), [1.into(), 2.into()]);
//     let g3: SparseTensor<Complex<f64>> = ufo::gamma(2.into(), [2.into(), 0.into()]);

//     let c = g1.contract(&g2).unwrap().contract(&g3).unwrap();
//     assert_eq!(
//         Vec::<Complex<f64>>::new(),
//         c.data(),
//         "Odd traces must vanish"
//     );

//     let d: SparseTensor<Complex<f64>> =
//         ufo::gamma(0.into(), [0.into(), 0.into()]).internal_contract();

//     assert_eq!(Vec::<Complex<f64>>::new(), d.data(), "Gammas are traceless");
// }

#[test]
fn matches() {
    let structur_a: HistoryStructure<String, ()> =
        HistoryStructure::from(NamedStructure::from_iter(
            [
                LibraryRep::from(Lorentz {}).slot(2, 3),
                Lorentz {}.slot(3, 2).cast(),
                Euclidean {}.slot(2, 2).cast(),
                Lorentz {}.slot(2, 1).cast(),
            ],
            "a".into(),
            None,
        ));
    let structur_b = HistoryStructure::from(NamedStructure::from_iter(
        [
            LibraryRep::from(Lorentz {}).dual().slot(2, 1),
            Lorentz {}.dual().slot(2, 3).cast(),
            Lorentz {}.dual().slot(2, 2).cast(),
            Euclidean {}.slot(2, 1).cast(),
        ],
        "b".into(),
        None,
    ));

    let a = structur_a.match_index(&structur_b);

    assert_eq!(Some((false, 3, 0)), a);
}

#[test]
fn mixed_tensor_contraction() {
    let im = Complex::new(1.5, 1.25);
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let structur_a: HistoryStructure<String, ()> =
        HistoryStructure::from(NamedStructure::from_iter(
            [
                LibraryRep::from(Euclidean {}).slot(2, 2),
                Euclidean {}.slot(2, 1).cast(),
            ],
            "a".into(),
            None,
        ));

    let a = SparseTensor::from_data(data_a.clone(), structur_a.clone()).unwrap();

    let structur_b = HistoryStructure::from(NamedStructure::from_iter(
        [
            LibraryRep::from(Euclidean {}).slot(2, 2),
            Euclidean {}.slot(2, 4).cast(),
        ],
        "b".into(),
        None,
    ));

    let b = DenseTensor::from_data(
        vec![
            im.mul_fallible(&1.0).unwrap(),
            2.0.mul_fallible(&im).unwrap(),
            3.0.mul_fallible(&im).unwrap(),
            4.0.mul_fallible(&im).unwrap(),
        ],
        structur_b.clone(),
    )
    .unwrap();

    let f = b.contract(&a).unwrap();

    assert_eq!(
        f.data,
        [
            1.0.mul_fallible(&im).unwrap(),
            6.0.mul_fallible(&im).unwrap(),
            2.0.mul_fallible(&im).unwrap(),
            8.0.mul_fallible(&im).unwrap()
        ]
    );

    let data_a = [
        (vec![0, 0], 1.0.mul_fallible(&im).unwrap()),
        (vec![1, 1], 2.0.mul_fallible(&im).unwrap()),
    ];

    let a = SparseTensor::from_data(data_a.clone(), structur_a).unwrap();

    let b = DenseTensor::from_data(vec![1.0, 2.0, 3.0, 4.0], structur_b).unwrap();

    let f = a.contract(&b).unwrap();
    assert_eq!(
        f.data,
        [
            1.0.mul_fallible(&im).unwrap(),
            2.0.mul_fallible(&im).unwrap(),
            6.0.mul_fallible(&im).unwrap(),
            8.0.mul_fallible(&im).unwrap()
        ]
    );
}

// #[test]
// fn tensor_net() {
//     let a: RealOrComplexTensor<f64, VecStructure> =
//         RealOrComplexTensor::Complex(DataTensor::from(ufo::gamma(1.into(), [2.into(), 3.into()])));
//     let b: RealOrComplexTensor<f64, _> =
//         RealOrComplexTensor::Complex(DataTensor::from(ufo::gamma(2.into(), [3.into(), 4.into()])));
//     let c: RealOrComplexTensor<f64, _> =
//         RealOrComplexTensor::Complex(DataTensor::from(ufo::gamma(3.into(), [4.into(), 5.into()])));
//     let d: RealOrComplexTensor<f64, _> =
//         RealOrComplexTensor::Complex(DataTensor::from(ufo::gamma(4.into(), [5.into(), 2.into()])));
//     let p: RealOrComplexTensor<f64, _> = RealOrComplexTensor::Real(DataTensor::from(
//         mink_four_vector(AbstractIndex::from(-1), [2., 3., 2., 1.]),
//     ));
//     let q: RealOrComplexTensor<f64, _> = RealOrComplexTensor::Real(DataTensor::from(
//         mink_four_vector(AbstractIndex::from(-2), [2., 3., 2., 1.]),
//     ));
//     let r: RealOrComplexTensor<f64, _> = RealOrComplexTensor::Real(DataTensor::from(
//         mink_four_vector(AbstractIndex::from(-3), [2., 3., 2., 1.]),
//     ));
//     let s: RealOrComplexTensor<f64, _> = RealOrComplexTensor::Real(DataTensor::from(
//         mink_four_vector(AbstractIndex::from(-4), [2., 3., 2., 1.]),
//     ));

//     let pslash = a.contract(&p).unwrap();
//     let qslash = b.contract(&q).unwrap();
//     let rslash = c.contract(&r).unwrap();
//     let sslash = d.contract(&s).unwrap();

//     let e = pslash
//         .contract(&qslash)
//         .unwrap()
//         .contract(&rslash)
//         .unwrap()
//         .contract(&sslash)
//         .unwrap();
//     // .contract(&s)
//     // .unwrap();

//     assert_eq!(
//         Complex::new(400., 0.),
//         e.scalar().unwrap().try_into_complex().unwrap()
//     ); //.scalar().unwrap());

//     let mut n: TensorNetwork<RealOrComplexTensor<f64, _>, Complex<f64>> =
//         TensorNetwork::from(vec![a, b, c, p, q, d, r, s]);

//     assert!(n.graph.validate_neighbors());

//     // println!("{}", n.dot());

//     assert_eq!(16, n.graph.neighbors.len());

//     n.contract().unwrap();

//     assert_eq!(0, n.graph.neighbors.len());
//     assert_eq!(
//         Complex::new(400., 0.),
//         n.result().unwrap().0.try_as_complex().unwrap().data()[0]
//     )
// }

#[test]
fn matchslot() {
    let a: LibrarySlot = Lorentz {}.dual().slot(4, 2).cast();
    let b: LibrarySlot = Lorentz {}.slot(4, 2).cast();
    assert!(a.matches(&b))
}
#[test]
fn contract_spensor() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];
    let structur_a: HistoryStructure<String, ()> =
        HistoryStructure::from(NamedStructure::from_iter(
            [
                LibraryRep::slot(Euclidean {}.into(), 2, 2),
                Euclidean {}.slot(2, 1).cast(),
            ],
            "a".into(),
            None,
        ));

    let a = SparseTensor::from_data(data_a.clone(), structur_a).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_b = HistoryStructure::from(NamedStructure::from_iter(
        [
            LibraryRep::slot(Euclidean {}.into(), 2, 1),
            Euclidean {}.slot(2, 3).cast(),
        ],
        "b".into(),
        None,
    ));

    let b = SparseTensor::from_data(data_b.clone(), structur_b).unwrap();

    let f = a.contract(&b).unwrap();

    let result = IndexMap::from([(vec![0, 1].into(), 2.0), (vec![1, 0].into(), 2.0)]);

    assert_eq!(f.hashmap(), result)
}

#[test]
fn sparse_addition() {
    let data_a = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_a: HistoryStructure<String, ()> =
        HistoryStructure::from(NamedStructure::from_iter(
            [
                LibraryRep::slot(Euclidean {}.into(), 2, 2),
                LibraryRep::slot(Euclidean {}.into(), 2, 1),
            ],
            "a".into(),
            None,
        ));

    let a = SparseTensor::from_data(data_a.clone(), structur_a).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_b = HistoryStructure::from(NamedStructure::from_iter(
        [
            LibraryRep::slot(Euclidean {}.into(), 2, 1),
            LibraryRep::slot(Euclidean {}.into(), 2, 2),
        ],
        "b".into(),
        None,
    ));

    let b = SparseTensor::from_data(data_b.clone(), structur_b).unwrap();

    let f = a.add_fallible(&b).unwrap();

    let result = IndexMap::from([(vec![0, 1].into(), 3.0), (vec![1, 0].into(), 3.0)]);

    assert_eq!(f.hashmap(), result)
}

#[test]
fn sparse_sub() {
    let data_a = [(vec![1, 0], 1.0), (vec![0, 1], 2.0)];
    let structur_a: HistoryStructure<String, ()> =
        HistoryStructure::from(NamedStructure::from_iter(
            [
                LibraryRep::slot(Euclidean {}.into(), 2, 2),
                LibraryRep::slot(Euclidean {}.into(), 2, 1),
            ],
            "a".into(),
            None,
        ));

    let a = SparseTensor::from_data(data_a.clone(), structur_a).unwrap();

    let data_b = [(vec![1, 0], 1.0), (vec![0, 1], 3.0)];

    let structur_b = HistoryStructure::from(NamedStructure::from_iter(
        [
            LibraryRep::slot(Euclidean {}.into(), 2, 2),
            LibraryRep::slot(Euclidean {}.into(), 2, 1),
        ],
        "a".into(),
        None,
    ));

    let b = SparseTensor::from_data(data_b.clone(), structur_b).unwrap();

    let f = a.sub_fallible(&b).unwrap();

    let result = IndexMap::from([(vec![0, 1].into(), -1.0)]);
    assert_eq!(f.hashmap(), result);
    // println!("{:?}", f);
}

#[test]
fn arithmetic_data() {
    let sa = test_structure(3, 1);
    let range = Some((-1000, 1000));

    let a: DataTensor<i32> = test_tensor(sa.clone(), 1, range).into();

    let b: DataTensor<i32> = test_tensor(sa.clone(), 2, range).into();

    println!("{:?}", a);
    println!("{:?}", b);

    let c = a.add_fallible(&b).unwrap();

    assert_ron_snapshot!(c.to_bare_dense().data());

    // let syma = sa.clone().shadow_with("a".into_id());
    // let symb = sa.clone().shadow_with("b".into_id());

    //  let symc = &syma + &symb;
}

#[test]
fn contract_densor_with_spensor() {
    let data_a = [(vec![0, 0], 1.0), (vec![1, 1], 2.0)];

    let structur_a: HistoryStructure<String, ()> =
        HistoryStructure::from(NamedStructure::from_iter(
            [
                LibraryRep::slot(Euclidean {}.into(), 2, 2),
                LibraryRep::slot(Euclidean {}.into(), 2, 1),
            ],
            "a".into(),
            None,
        ));

    let a = SparseTensor::from_data(data_a.clone(), structur_a).unwrap();

    let data_b = [1.0, 2.0, 3.0, 4.0];
    let structur_b = HistoryStructure::from(NamedStructure::from_iter(
        [
            LibraryRep::slot(Euclidean {}.into(), 2, 1),
            LibraryRep::slot(Euclidean {}.into(), 2, 4),
        ],
        "b".into(),
        None,
    ));

    let b = DenseTensor::from_data(data_b.to_vec(), structur_b).unwrap();

    let f = a.contract(&b).unwrap();

    assert_eq!(f.data, [1.0, 2.0, 6.0, 8.0]);
}

// #[test]
// fn symbolic_zeros() {
//     let mut state = State::get_global_state().write().unwrap();
//     let ws = Workspace::new();
//     let structure = TensorSkeleton::from(NamedStructure::from_iter([(1, 2), (3, 2)]seTen,None)sor::symbolic_zeros(structure.clone());

//     let zeros: DenseTensor<f64> = DenseTensor::default(structure);

//     assert_eq!(sym_zeros, zeros.to_symbolic(&ws, &mut state));
// }

#[test]
#[cfg(feature = "shadowing")]
fn evaluate() {
    use crate::{parametric::atomcore::TensorAtomMaps, shadowing::Shadowable};

    let structure: NamedStructure<String, ()> = test_structure(3, 1).to_named("a".into(), None);

    let a: TensorShell<_> = structure.clone().into();

    let a = a.expanded_shadow().unwrap();

    let adata = test_tensor(structure, 1, Some((-100., 100.))).to_dense();

    let mut const_map: HashMap<AtomView<'_>, f64> = HashMap::new();

    a.append_const_map(&adata, &mut const_map);

    let fn_map = HashMap::new();

    let aev: DenseTensor<f64, _> = a.evaluate(|r| r.into(), &const_map, &fn_map).unwrap();

    assert_eq!(aev.data(), adata.data());
}

#[test]
#[cfg(feature = "shadowing")]
fn convert_sym() {
    let i: Complex<f64> = Complex::new(0.0, 1.0);
    let mut data_b: Vec<Complex<f64>> = vec![i * Complex::from(5.0), Complex::<f64>::from(2.6) + i];
    data_b.append(
        &mut [3.34, -17.125, 5.0, 6.0]
            .iter()
            .map(|x| Complex::from(*x))
            .collect::<Vec<_>>(),
    );
    let structur_b: HistoryStructure<String, ()> =
        HistoryStructure::from(NamedStructure::from_iter(
            [
                LibraryRep::slot(Euclidean {}.into(), 2, 1),
                LibraryRep::slot(Euclidean {}.into(), 3, 4),
            ],
            "b".into(),
            None,
        ));
    let b = DenseTensor::from_data(data_b.to_vec(), structur_b).unwrap();

    let symb: DenseTensor<Atom, _> = b.try_into_upgrade().unwrap();

    let expected_data: Vec<Atom> = [
        "5*𝑖",
        "𝑖+5854679515581645/2251799813685248",
        "940126422213591/281474976710656",
        "-137/8",
        "5",
        "6",
    ]
    .iter()
    .map(|x| parse!(x).unwrap())
    .collect();

    assert_eq!(
        symb.iter_expanded()
            .map(|(_, x)| x.clone())
            .collect::<Vec<_>>(),
        expected_data
    );
}

#[test]
#[cfg(feature = "shadowing")]
fn simple_multi_contract_sym() {
    let structa = VecStructure::from_iter([
        Euclidean {}.slot(3, 1),
        Euclidean {}.slot(4, 2),
        Euclidean {}.slot(4, 3),
    ]);
    // let structa = structa.to_named("a");
    let structb = VecStructure::from_iter([
        Euclidean {}.slot(4, 2),
        Euclidean {}.slot(3, 4),
        Euclidean {}.slot(4, 3),
    ]);
    // let structb = structb.to_named("b");

    let _a = DenseTensor::from_data(
        vec![
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 12, 13, 14, 15, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11,
            12, 13, 14, 15, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 12, 13, 14, 15, 16,
        ],
        structa.clone(),
    )
    .unwrap();

    let _b = DenseTensor::from_data(
        vec![
            3, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 3, 2, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0,
            3, 2, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 3, 2, 1,
        ],
        structb.clone(),
    )
    .unwrap();

    let a: DataTensor<Atom, NamedStructure<String, ()>> = structa
        .to_named("a".into(), None)
        .to_dense_expanded_labels()
        .unwrap()
        .into();
    let b: DataTensor<Atom, NamedStructure<String, ()>> = structb
        .to_named("b".into(), None)
        .to_dense_expanded_labels()
        .unwrap()
        .into();

    let _f = a.contract(&b).unwrap();
}

#[test]
fn empty_densor() {
    let empty_structure = Vec::<Slot<_>>::new();

    let empty: DenseTensor<f64> = DenseTensor::default(empty_structure.into());

    assert_eq!(*empty.get_ref(&[]).unwrap(), 0.0);
}

#[test]
fn complex() {
    let _structur = test_structure(2, 1);

    let _r = Complex::new(1.0, 2.0);
    let _p = Complex::new(3.0, 4.0);
}

#[test]
#[cfg(feature = "shadowing")]
fn symbolic_contract() {
    let structura: NamedStructure<String, ()> = NamedStructure::from_iter(
        [
            LibraryRep::slot(Euclidean {}.into(), 2, 1),
            LibraryRep::slot(Euclidean {}.into(), 3, 4),
        ],
        "T".into(),
        None,
    );

    let structurb: NamedStructure<String, ()> = NamedStructure::from_iter(
        [
            LibraryRep::slot(Euclidean {}.into(), 2, 3),
            LibraryRep::slot(Euclidean {}.into(), 3, 2),
        ],
        "P".into(),
        None,
    );

    let mink = Lorentz {}.rep(4);
    let mu: Slot<LibraryRep> = mink.slot(0).cast();
    let bis = Euclidean {}.rep(4);

    let i = bis.slot(1).cast();
    let j = bis.slot(2).cast();
    let k = bis.slot(9).cast();

    let _structure: NamedStructure<String, ()> =
        NamedStructure::from_iter([mu, i, j], "γ".into(), None);
    let _p_struct: NamedStructure<String, ()> = NamedStructure::from_iter([mu], "p".into(), None);
    let _t_struct: NamedStructure<String, ()> =
        NamedStructure::from_iter([i, j, k], "T".into(), None);

    let a = SymbolicTensor::from_named(&structura.to_shell()).unwrap();
    let b = SymbolicTensor::from_named(&structurb.to_shell()).unwrap();
    let f = a.contract(&b).unwrap();

    println!("{}", f);

    assert_eq!(
        *f.get_atom(),
        parse!("T(euc(2,1),euc(3,4))*P(euc(2,3),euc(3,2))").unwrap()
    );

    let a = f.to_network(&EXPLICIT_TENSOR_MAP.read().unwrap()).unwrap();

    // let syms = a.to_symbolic_tensor_vec();

    // for s in syms {
    //     println!("{:?}", s.structure());
    // }

    println!("{}", a.dot());
}

#[test]
fn test_fallible_mul() {
    let a: i32 = 4;
    let b: f64 = 4.;
    let mut c: f64 = a.mul_fallible(&b).unwrap();
    c.add_assign_fallible(&a);
    let d: Option<f64> = b.mul_fallible(&a);
    let a: &i32 = &a;
    let e: Option<f64> = a.mul_fallible(&b);
    assert_eq!(c, 20.);
    assert_eq!(d, Some(16.));
    assert_eq!(e, Some(16.));

    #[cfg(feature = "shadowing")]
    {
        let a = &parse!("a(2)").unwrap();

        let b = &parse!("b(1)").unwrap();
        let i = Atom::var(Atom::I);
        let mut f = a.mul_fallible(&4.).unwrap();
        f.add_assign_fallible(b);
        f.add_assign_fallible(&i);

        let function_map = HashMap::new();

        let mut const_map = HashMap::new();
        const_map.insert(i.as_view(), Complex::<f64>::new(0., 1.).into());

        const_map.insert(a.as_view(), Complex::<f64>::new(3., 1.).into());

        const_map.insert(b.as_view(), Complex::<f64>::new(3., 1.).into());

        let ev: symbolica::domains::float::Complex<f64> =
            f.evaluate(|r| r.into(), &const_map, &function_map).unwrap();

        println!("{}", ev);
        // print!("{}", f.unwrap());

        let g = Complex::new(0.1, 3.);

        let mut h = a.sub_fallible(&g).unwrap();

        h.add_assign_fallible(a);
        let _f = a.mul_fallible(a);

        println!("{}", h);
    }
}

// #[test]
// #[cfg(feature = "shadowing")]
// fn get_license_key() {
//     symbolica::LicenseManager::new();
// }

#[test]
fn duals() {
    let a = Lorentz {}.dual().rep(4);
    let mu = a.slot(1);
    let nu = a.dual().slot(1);

    assert!(mu.matches(&nu));

    let mu: LibrarySlot = mu.cast();

    assert!(mu.matches(&nu.cast()));

    let rho: LibrarySlot = Minkowski {}.dual().slot(4, 1).cast();

    assert!(!mu.matches(&rho))
}

#[test]
fn transpose() {
    let a = DenseTensor::from_data(
        vec![1, 2],
        VecStructure::<Euclidean>::from_iter([Euclidean {}.slot(2, 1)]),
    )
    .unwrap();

    let b = DenseTensor::from_data(
        vec![3, 4],
        VecStructure::from_iter([Euclidean {}.slot(2, 2)]),
    )
    .unwrap();
    // let m = DenseTensor::from_data(
    //     vec![0, 1, 0, 0],
    //     VecStructure::new(vec![
    //         Bispinor::new_slot_selfless(2, 2).into(),
    //         Bispinor::new_slot_selfless(2, 1).into(),
    //     ]),
    // )
    // .unwrap();

    let mtranspose = DenseTensor::from_data(
        vec![0, 1, 0, 0],
        VecStructure::from_iter([Euclidean {}.slot(2, 1), Euclidean {}.slot(2, 2)]),
    )
    .unwrap();

    let mdatatransposed = DenseTensor::from_data(
        vec![0, 0, 1, 0],
        VecStructure::new(vec![Euclidean {}.slot(2, 2), Euclidean {}.slot(2, 1)]),
    )
    .unwrap();

    let ftranspose = a.contract(&mtranspose).unwrap().contract(&b).unwrap();

    let fdatatransposed = a.contract(&mdatatransposed).unwrap().contract(&b).unwrap();

    assert_eq!(fdatatransposed.data(), ftranspose.data());
}
