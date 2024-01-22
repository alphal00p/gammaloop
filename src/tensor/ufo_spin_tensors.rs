use super::{
    Fouroot, Representation, Representation::Euclidean, Representation::Lorentz, SparseTensor,
    TensorStructure, VecSlotExtension,
};
use num::{Complex, Float, One, Zero};

#[allow(dead_code)]
pub fn identity<T>(indices: (usize, usize), signature: Representation) -> SparseTensor<Complex<T>>
where
    T: One + Zero,
{
    //TODO: make it just swap indices
    let structure = TensorStructure::from_idxsing(&[indices.0, indices.1], &[signature, signature]);
    let mut identity = SparseTensor::empty(structure);
    for i in 0..signature.into() {
        identity
            .set(&[i, i], Complex::<T>::new(T::one(), T::zero()))
            .unwrap();
    }
    identity
}

#[allow(dead_code)]
pub fn lorentz_identity<T>(indices: (usize, usize)) -> SparseTensor<Complex<T>>
where
    T: One + Zero,
{
    // IdentityL(1,2) (Lorentz) Kronecker delta δ^μ1_μ1
    let signature = Lorentz(4);
    identity(indices, signature)
}

#[allow(dead_code)]
pub fn euclidean_identity<T>(indices: (usize, usize)) -> SparseTensor<Complex<T>>
where
    T: One + Zero,
{
    // Identity(1,2) (Spinorial) Kronecker delta δ_s1_s2
    let signature = Euclidean(4);
    identity(indices, signature)
}

#[allow(dead_code)]
pub fn gamma(minkindex: usize, indices: (usize, usize)) -> SparseTensor<Complex<i8>> {
    // Gamma(1,2,3) Dirac matrix (γ^μ1)_s2_s3
    let structure = TensorStructure::from_idxsing(
        &[indices.0, indices.1, minkindex],
        &[Euclidean(4), Euclidean(4), Lorentz(4)],
    );

    let c1 = Complex::<i8>::new(1, 0);
    let cn1 = Complex::<i8>::new(-1, 0);
    let ci = Complex::<i8>::new(0, 1);
    let cni = Complex::<i8>::new(0, -1);

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

pub fn gamma5<T>(indices: (usize, usize)) -> SparseTensor<Complex<T>>
where
    T: One + Zero + Copy,
{
    let structure =
        TensorStructure::from_idxsing(&[indices.0, indices.1], &[Euclidean(4), Euclidean(4)]);

    let c1 = Complex::<T>::new(T::one(), T::zero());

    let mut gamma5 = SparseTensor::empty(structure);

    gamma5.set(&[0, 2], c1).unwrap();
    gamma5.set(&[1, 3], c1).unwrap();
    gamma5.set(&[2, 0], c1).unwrap();
    gamma5.set(&[3, 1], c1).unwrap();

    gamma5
}

pub fn proj_m<T>(indices: (usize, usize)) -> SparseTensor<Complex<T>>
where
    T: Float,
{
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
    let structure =
        TensorStructure::from_idxsing(&[indices.0, indices.1], &[Euclidean(4), Euclidean(4)]);

    let chalf = Complex::<T>::new(T::from(0.5).unwrap(), T::zero());
    let cnhalf = Complex::<T>::new(T::from(-0.5).unwrap(), T::zero());

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

pub fn proj_p<T>(indices: (usize, usize)) -> SparseTensor<Complex<T>>
where
    T: Float,
{
    // ProjP(1,2) Right chirality projector (( 1+γ5)/ 2 )_s1_s2
    let structure =
        TensorStructure::from_idxsing(&[indices.0, indices.1], &[Euclidean(4), Euclidean(4)]);

    let chalf = Complex::<T>::new(T::from(0.5).unwrap(), T::zero());

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

pub fn sigma<T>(indices: (usize, usize), minkdices: (usize, usize)) -> SparseTensor<Complex<T>>
where
    T: One + Zero + std::ops::Neg<Output = T> + Copy,
{
    let structure = TensorStructure::from_idxsing(
        &[indices.0, indices.1, minkdices.0, minkdices.1],
        &[Euclidean(4), Euclidean(4), Lorentz(4), Lorentz(4)],
    );

    let c1 = Complex::<T>::new(T::one(), T::zero());
    let cn1 = Complex::<T>::new(-T::one(), T::zero());
    let ci = Complex::<T>::new(T::zero(), T::one());
    let cni = Complex::<T>::new(T::zero(), -T::one());

    let mut sigma = SparseTensor::empty(structure);
    sigma.set(&[0, 2, 0, 1], c1).unwrap();
    sigma.set(&[0, 2, 3, 0], c1).unwrap();
    sigma.set(&[0, 3, 1, 2], c1).unwrap();
    sigma.set(&[1, 0, 2, 2], c1).unwrap();
    sigma.set(&[1, 1, 1, 2], c1).unwrap();
    sigma.set(&[1, 3, 0, 2], c1).unwrap();
    sigma.set(&[2, 2, 1, 0], c1).unwrap();
    sigma.set(&[2, 2, 2, 1], c1).unwrap();
    sigma.set(&[2, 3, 3, 2], c1).unwrap();
    sigma.set(&[3, 0, 0, 2], c1).unwrap();
    sigma.set(&[3, 3, 2, 2], c1).unwrap();
    sigma.set(&[3, 1, 3, 2], c1).unwrap();
    sigma.set(&[0, 1, 3, 0], ci).unwrap();
    sigma.set(&[0, 3, 1, 1], ci).unwrap();
    sigma.set(&[0, 3, 2, 0], ci).unwrap();
    sigma.set(&[1, 0, 3, 3], ci).unwrap();
    sigma.set(&[1, 1, 0, 3], ci).unwrap();
    sigma.set(&[1, 1, 2, 0], ci).unwrap();
    sigma.set(&[2, 1, 1, 0], ci).unwrap();
    sigma.set(&[2, 3, 0, 0], ci).unwrap();
    sigma.set(&[2, 3, 3, 1], ci).unwrap();
    sigma.set(&[3, 0, 1, 3], ci).unwrap();
    sigma.set(&[3, 1, 0, 0], ci).unwrap();
    sigma.set(&[3, 1, 2, 3], ci).unwrap();
    sigma.set(&[0, 0, 3, 2], cn1).unwrap();
    sigma.set(&[0, 1, 0, 2], cn1).unwrap();
    sigma.set(&[0, 2, 1, 3], cn1).unwrap();
    sigma.set(&[1, 2, 0, 3], cn1).unwrap();
    sigma.set(&[1, 2, 1, 1], cn1).unwrap();
    sigma.set(&[1, 2, 2, 0], cn1).unwrap();
    sigma.set(&[2, 0, 1, 2], cn1).unwrap();
    sigma.set(&[2, 1, 2, 2], cn1).unwrap();
    sigma.set(&[2, 2, 3, 3], cn1).unwrap();
    sigma.set(&[3, 2, 0, 0], cn1).unwrap();
    sigma.set(&[3, 2, 2, 3], cn1).unwrap();
    sigma.set(&[3, 2, 3, 1], cn1).unwrap();
    sigma.set(&[0, 0, 2, 3], cni).unwrap();
    sigma.set(&[0, 0, 3, 1], cni).unwrap();
    sigma.set(&[0, 1, 1, 3], cni).unwrap();
    sigma.set(&[1, 0, 2, 1], cni).unwrap();
    sigma.set(&[1, 3, 0, 1], cni).unwrap();
    sigma.set(&[1, 3, 3, 0], cni).unwrap();
    sigma.set(&[2, 0, 0, 3], cni).unwrap();
    sigma.set(&[2, 0, 1, 1], cni).unwrap();
    sigma.set(&[2, 1, 3, 3], cni).unwrap();
    sigma.set(&[3, 0, 0, 1], cni).unwrap();
    sigma.set(&[3, 3, 1, 0], cni).unwrap();
    sigma.set(&[3, 3, 2, 1], cni).unwrap();

    sigma
}
