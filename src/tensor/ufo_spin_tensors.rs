use super::{
    ContractableWithSparse, Representation, Representation::Euclidean, Representation::Lorentz,
    SparseTensor, TensorStructure, VecSlotExtension,
};
use num::complex::Complex64;

#[allow(dead_code)]
pub fn identity(indices: (usize, usize), signature: Representation) -> SparseTensor<Complex64> {
    //TODO: make it just swap indices
    let structure = TensorStructure::from_idxsing(&[indices.0, indices.1], &[signature, signature]);
    let mut identity = SparseTensor::empty(structure);
    for i in 0..signature.into() {
        identity.set(&[i, i], Complex64::new(1.0, 0.0)).unwrap();
    }
    identity
}

#[allow(dead_code)]
pub fn lorentz_identity(indices: (usize, usize)) -> SparseTensor<Complex64> {
    // IdentityL(1,2) (Lorentz) Kronecker delta δ^μ1_μ1
    let signature = Lorentz(4);
    identity(indices, signature)
}

#[allow(dead_code)]
pub fn euclidean_identity(indices: (usize, usize)) -> SparseTensor<Complex64> {
    // Identity(1,2) (Spinorial) Kronecker delta δ_s1_s2
    let signature = Euclidean(4);
    identity(indices, signature)
}

#[allow(dead_code)]
pub fn gamma(minkindex: usize, indices: (usize, usize)) -> SparseTensor<Complex64> {
    // Gamma(1,2,3) Dirac matrix (γ^μ1)_s2_s3
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

pub fn gamma5(indices: (usize, usize)) -> SparseTensor<Complex64> {
    let structure =
        TensorStructure::from_idxsing(&[indices.0, indices.1], &[Euclidean(4), Euclidean(4)]);

    let c1 = Complex64::new(1.0, 0.0);

    let mut gamma5 = SparseTensor::empty(structure);

    gamma5.set(&[0, 2], c1).unwrap();
    gamma5.set(&[1, 3], c1).unwrap();
    gamma5.set(&[2, 0], c1).unwrap();
    gamma5.set(&[3, 1], c1).unwrap();

    gamma5
}

pub fn proj_m(indices: (usize, usize)) -> SparseTensor<Complex64> {
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
    let structure =
        TensorStructure::from_idxsing(&[indices.0, indices.1], &[Euclidean(4), Euclidean(4)]);

    let chalf = Complex64::new(0.5, 0.0);
    let cnhalf = Complex64::new(-0.5, 0.0);

    let mut proj_m = SparseTensor::empty(structure);

    proj_m.set(&[0, 0], chalf).unwrap();
    proj_m.set(&[1, 1], chalf).unwrap();
    proj_m.set(&[2, 2], chalf).unwrap();
    proj_m.set(&[3, 3], chalf).unwrap();

    proj_m.set(&[0, 2], cnhalf).unwrap();
    proj_m.set(&[1, 3], cnhalf).unwrap();
    proj_m.set(&[2, 0], cnhalf).unwrap();
    proj_m.set(&[3, 1], cnhalf).unwrap();

    proj_m
}

pub fn proj_p(indices: (usize, usize)) -> SparseTensor<Complex64> {
    // ProjP(1,2) Right chirality projector (( 1+γ5)/ 2 )_s1_s2
    let structure =
        TensorStructure::from_idxsing(&[indices.0, indices.1], &[Euclidean(4), Euclidean(4)]);

    let chalf = Complex64::new(0.5, 0.0);

    let mut proj_p = SparseTensor::empty(structure);

    proj_p.set(&[0, 0], chalf).unwrap();
    proj_p.set(&[1, 1], chalf).unwrap();
    proj_p.set(&[2, 2], chalf).unwrap();
    proj_p.set(&[3, 3], chalf).unwrap();

    proj_p.set(&[0, 2], chalf).unwrap();
    proj_p.set(&[1, 3], chalf).unwrap();
    proj_p.set(&[2, 0], chalf).unwrap();
    proj_p.set(&[3, 1], chalf).unwrap();

    proj_p
}

pub fn sigma(indices: (usize, usize), minkdices: (usize, usize)) -> SparseTensor<Complex64> {
    let k = indices.0 + indices.1;

    let imhalf = Complex64::new(0.0, 0.5);

    let gamma1 = gamma(minkdices.0, (indices.0, k));
    let gamma2 = gamma(minkdices.1, (k, indices.1));
    let gamma3 = gamma(minkdices.1, (indices.0, k));
    let gamma4 = gamma(minkdices.0, (k, indices.1));

    (gamma1.contract_with_sparse(&gamma2).unwrap() - gamma3.contract_with_sparse(&gamma4).unwrap())
        * imhalf
}
