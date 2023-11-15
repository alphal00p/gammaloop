// Gamma chain example

use std::{
    ops::Neg,time::Instant,
};

use _gammaloop::tensor::{
    DenseTensor,
    Signature::Lorentz,
    Signature::{self, Euclidean},
    SparseTensor, TensorStructure, VecSlotExt,
};

use num_complex::Complex64;
use num_traits::{Num, ToPrimitive};

fn gamma(minkindex: usize, indices: (usize, usize)) -> SparseTensor<Complex64> {
    let structure = TensorStructure::from_idxsing(
        &[indices.0, indices.1, minkindex],
        &[Euclidean(4), Euclidean(4), Lorentz(4)],
    );

    let c1 = Complex64::new(1.0, 0.0);
    let cn1 = Complex64::new(-1.0, 0.0);
    let ci = Complex64::new(0.0, 1.0);
    let cni = Complex64::new(0.0, -1.0);

    let mut gamma = SparseTensor::empty(structure);

    // dirac gamma matrices

    gamma.set(&[0, 0, 0], c1).unwrap();
    gamma.set(&[1, 1, 0], c1).unwrap();
    gamma.set(&[2, 2, 0], cn1).unwrap();
    gamma.set(&[3, 3, 0], cn1).unwrap();

    gamma.set(&[0, 3, 1], c1).unwrap();
    gamma.set(&[1, 2, 1], c1).unwrap();
    gamma.set(&[2, 1, 1], cn1).unwrap();
    gamma.set(&[3, 0, 1], cn1).unwrap();

    gamma.set(&[0, 3, 2], cni).unwrap();
    gamma.set(&[1, 2, 2], ci).unwrap();
    gamma.set(&[2, 1, 2], ci).unwrap();
    gamma.set(&[3, 0, 2], cni).unwrap();

    gamma.set(&[0, 2, 3], c1).unwrap();
    gamma.set(&[1, 3, 3], cn1).unwrap();
    gamma.set(&[2, 0, 3], cn1).unwrap();
    gamma.set(&[3, 1, 3], c1).unwrap();

    gamma
}

fn pslash(indices: (usize, usize), p: [Complex64; 4]) -> DenseTensor<Complex64> {


    let minkindex = indices.0 + indices.1;

    let p = DenseTensor::from_data(
        &p,
        TensorStructure::from_idxsing(&[minkindex], &[Lorentz(4)]),
    )
    .unwrap();

    let  pslash = gamma(minkindex, indices).contract_with_dense(&p).unwrap();

    pslash
}


#[allow(dead_code)]
fn mink_four_vector<T>(index: usize, p: [T; 4]) -> DenseTensor<T>
where
    T: Num + std::default::Default + std::clone::Clone,
{
    DenseTensor::from_data(&p, TensorStructure::from_idxsing(&[index], &[Lorentz(4)])).unwrap()
}

fn eucl_four_vector<T>(index: usize, p: [T; 4]) -> DenseTensor<T>
where
    T: Num + std::default::Default + std::clone::Clone,
{
 DenseTensor::from_data(&p, TensorStructure::from_idxsing(&[index], &[Euclidean(4)])).unwrap()
}

fn benchmark_chain(
    minkindices: &[i32],
    vbar: [Complex64; 4],
    u: [Complex64; 4],
) -> DenseTensor<Complex64> {
    let mut i = 0;
    let mut contracting_index = 0;
    let mut result = eucl_four_vector(contracting_index, vbar);
    for m in minkindices {
        if *m > 0 {
            let p = [
                Complex64::new(1.0 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.1 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.2 + 0.01 * i.to_f64().unwrap(), 0.0),
                Complex64::new(1.3 + 0.01 * i.to_f64().unwrap(), 0.0),
            ];
            i += 1;
            let pslash = pslash((contracting_index, contracting_index + 1), p);
            result = pslash.contract_with_dense(&result).unwrap();
        } else {
            result = gamma(
                usize::try_from(m.neg()).unwrap(),
                (contracting_index, contracting_index + 1),
            )
            .contract_with_dense(&result)
            .unwrap();
        }
        contracting_index += 1;
    }
    result
        .contract_with_dense(&eucl_four_vector(contracting_index, u))
        .unwrap()
}


#[allow(dead_code)]
fn identity(indices: (usize, usize), signature: Signature) -> DenseTensor<Complex64> {
    let structure = TensorStructure::from_idxsing(&[indices.0, indices.1], &[signature, signature]);
    let mut identity = DenseTensor::default(structure);
    for i in 0..signature.into() {
        identity.set(&[i, i], Complex64::new(1.0, 0.0));
    }
    identity
}

#[allow(unused_variables)]
fn main() {
    // let p = [Complex64::new(1.0, 0.0); 4];
    // let p1 = pslash((1, 2), p);

    // let trg = gamma(1, (2, 2)).internal_contract();

    let one = Complex64::new(1.0, 0.0);
    let zero = Complex64::new(0.0, 0.0);

    // let gammamu = gamma(1, (2, 3))
    //     .contract_with_dense(&identity((3, 4), Euclidean(4)))
    //     .unwrap();

    // // println!("{:?}", gammamu);

    // let gammunu = gamma(1, (4, 5))
    //     .contract_with_dense(&identity((5, 6), Euclidean(4)))
    //     .unwrap();

    // println!("{:?}", gammamu.contract_with_dense(&gammunu).unwrap());

    let vbar = [one * 3.0, one * 3.1, one * 3.2, one * 3.3];
    let u = [one * 4.0, one * 4.1, one * 4.2, one * 4.3];
    // let a = eucl_four_vector(1, vbar);
    // let b = eucl_four_vector(1, u);
    // println!("{:?}", a.contract_with_dense(&b));

    // let p11 = pslash((1, 2), [one, zero, zero, zero])
    //     .contract_with_dense(&pslash((2, 1), [one, zero, zero, zero]));

    // println!("P {:?}", p11);

    let spacings: [i32; 2] = [30, 400];
    let mut start = 1;
    let mut ranges = Vec::new();

    for &spacing in spacings.iter() {
        ranges.push((start..start + spacing).chain(std::iter::once(-1)));
        start += spacing;
    }

    let vec: Vec<i32> = ranges.into_iter().flat_map(|x| x).collect();

    let start = Instant::now();
    let chain = benchmark_chain(&vec, vbar, u);
    let duration = start.elapsed();

    println!("{:?} in {:?}", chain, duration);
}
