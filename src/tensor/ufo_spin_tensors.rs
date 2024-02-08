use super::{
    DenseTensor, IntoId, Representation, Representation::Euclidean, Representation::Lorentz,
    SparseTensor, TensorSkeleton,
};
use num::{Complex, Float, One, Zero};


use symbolica::{
    representations::{Atom, Identifier},
    state::Workspace,
};

#[allow(dead_code)]
pub fn identity<T>(indices: (usize, usize), signature: Representation) -> SparseTensor<Complex<T>>
where
    T: One + Zero,
{
    //TODO: make it just swap indices
    let structure = TensorSkeleton::from_idxsing(
        &[(indices.0, signature), (indices.1, signature)],
        "id".into(),
    );
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

pub fn mink_four_vector<T>(index: usize, p: &[T; 4]) -> DenseTensor<T>
where
    T: Clone,
{
    DenseTensor::from_data(
        p,
        TensorSkeleton::from_idxsing(&[(index, Lorentz(4))], "p".into()),
    )
    .unwrap()
}

pub fn mink_four_vector_sym<T>(
    index: usize,
    p: &[T; 4],
    state: &mut symbolica::state::State,
) -> DenseTensor<T, Identifier>
where
    T: Clone,
{
    DenseTensor::from_data(
        p,
        TensorSkeleton::from_idxsing(
            &[(index, Lorentz(4))],
            state.get_or_insert_fn("p", None).unwrap(),
        ),
    )
    .unwrap()
}

pub fn euclidean_four_vector<T>(index: usize, p: &[T; 4]) -> DenseTensor<T>
where
    T: Clone,
{
    DenseTensor::from_data(
        p,
        TensorSkeleton::from_idxsing(&[(index, Euclidean(4))], "p".into()),
    )
    .unwrap()
}

pub fn euclidean_four_vector_sym<T>(
    index: usize,
    p: &[T; 4],
    state: &mut symbolica::state::State,
) -> DenseTensor<T, Identifier>
where
    T: Clone,
{
    DenseTensor::from_data(
        p,
        TensorSkeleton::from_idxsing(
            &[(index, Euclidean(4))],
            state.get_or_insert_fn("p", None).unwrap(),
        ),
    )
    .unwrap()
}

pub fn param_mink_four_vector<N>(
    index: usize,
    name: N,
    state: &mut symbolica::state::State,
    ws: &Workspace,
) -> DenseTensor<Atom, N>
where
    N: Clone + IntoId,
{
    TensorSkeleton::from_idxsing(&[(index, Lorentz(4))], name)
        .to_dense(state, ws)
        .unwrap()
}

pub fn param_euclidean_four_vector<N>(
    index: usize,
    name: N,
    state: &mut symbolica::state::State,
    ws: &Workspace,
) -> DenseTensor<Atom, N>
where
    N: Clone + IntoId,
{
    TensorSkeleton::from_idxsing(&[(index, Euclidean(4))], name)
        .to_dense(state, ws)
        .unwrap()
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
pub fn gamma<T>(minkindex: usize, indices: (usize, usize)) -> SparseTensor<Complex<T>>
where
    T: One + Zero + Copy + std::ops::Neg<Output = T>,
{
    // Gamma(1,2,3) Dirac matrix (γ^μ1)_s2_s3
    let structure = TensorSkeleton::from_idxsing(
        &[
            (indices.0, Euclidean(4)),
            (indices.1, Euclidean(4)),
            (minkindex, Lorentz(4)),
        ],
        "γ".into(),
    );

    gamma_data(structure)
}

pub fn gammasym<T>(
    minkindex: usize,
    indices: (usize, usize),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, Identifier>
where
    T: One + Zero + Copy + std::ops::Neg<Output = T>,
{
    let structure = TensorSkeleton::from_idxsing(
        &[
            (indices.0, Euclidean(4)),
            (indices.1, Euclidean(4)),
            (minkindex, Lorentz(4)),
        ],
        state.get_or_insert_fn("γ", None).unwrap(),
    );

    gamma_data(structure)
}

fn gamma_data<N, T>(structure: TensorSkeleton<N>) -> SparseTensor<Complex<T>, N>
where
    T: One + Zero + Copy + std::ops::Neg<Output = T>,
{
    let c1 = Complex::<T>::new(T::one(), T::zero());
    let cn1 = Complex::<T>::new(-T::one(), T::zero());
    let ci = Complex::<T>::new(T::zero(), T::one());
    let cni = Complex::<T>::new(T::zero(), -T::one());
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

    gamma //.to_dense()
}

pub fn gamma5<T>(indices: (usize, usize)) -> SparseTensor<Complex<T>>
where
    T: One + Zero + Copy,
{
    let structure = TensorSkeleton::from_idxsing(
        &[(indices.0, Euclidean(4)), (indices.1, Euclidean(4))],
        "γ5".into(),
    );

    gamma5_data(structure)
}

pub fn gamma5sym<T>(
    indices: (usize, usize),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, Identifier>
where
    T: One + Zero + Copy,
{
    let structure = TensorSkeleton::from_idxsing(
        &[
            (indices.0, Euclidean(4)),
            (indices.1, Euclidean(4)),
            (4, Lorentz(4)),
        ],
        state.get_or_insert_fn("γ5", None).unwrap(),
    );

    gamma5_data(structure)
}

fn gamma5_data<T, N>(structure: TensorSkeleton<N>) -> SparseTensor<Complex<T>, N>
where
    T: One + Zero + Copy,
{
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
    let structure = TensorSkeleton::from_idxsing(
        &[(indices.0, Euclidean(4)), (indices.1, Euclidean(4))],
        "ProjM".into(),
    );

    proj_m_data(structure)
}

pub fn proj_msym<T>(
    indices: (usize, usize),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, Identifier>
where
    T: Float,
{
    let structure = TensorSkeleton::from_idxsing(
        &[
            (indices.0, Euclidean(4)),
            (indices.1, Euclidean(4)),
            (4, Lorentz(4)),
        ],
        state.get_or_insert_fn("ProjM", None).unwrap(),
    );

    proj_m_data(structure)
}

fn proj_m_data<T, N>(structure: TensorSkeleton<N>) -> SparseTensor<Complex<T>, N>
where
    T: Float,
{
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
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
    let structure = TensorSkeleton::from_idxsing(
        &[(indices.0, Euclidean(4)), (indices.1, Euclidean(4))],
        "ProjP".into(),
    );

    proj_p_data(structure)
}

pub fn proj_psym<T>(
    indices: (usize, usize),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, Identifier>
where
    T: Float,
{
    let structure = TensorSkeleton::from_idxsing(
        &[
            (indices.0, Euclidean(4)),
            (indices.1, Euclidean(4)),
            (4, Lorentz(4)),
        ],
        state.get_or_insert_fn("ProjP", None).unwrap(),
    );

    proj_p_data(structure)
}

pub fn proj_p_data<T, N>(structure: TensorSkeleton<N>) -> SparseTensor<Complex<T>, N>
where
    T: Float,
{
    // ProjP(1,2) Right chirality projector (( 1+γ5)/ 2 )_s1_s2
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
    let structure = TensorSkeleton::from_idxsing(
        &[
            (indices.0, Euclidean(4)),
            (indices.1, Euclidean(4)),
            (minkdices.0, Lorentz(4)),
            (minkdices.1, Lorentz(4)),
        ],
        "σ".into(),
    );

    sigma_data(structure)
}

pub fn sigmasym<T>(
    indices: (usize, usize),
    minkdices: (usize, usize),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, Identifier>
where
    T: One + Zero + std::ops::Neg<Output = T> + Copy,
{
    let structure = TensorSkeleton::from_idxsing(
        &[
            (indices.0, Euclidean(4)),
            (indices.1, Euclidean(4)),
            (minkdices.0, Lorentz(4)),
            (minkdices.1, Lorentz(4)),
        ],
        state.get_or_insert_fn("σ", None).unwrap(),
    );

    sigma_data(structure)
}

fn sigma_data<T, N>(structure: TensorSkeleton<N>) -> SparseTensor<Complex<T>, N>
where
    T: One + Zero + std::ops::Neg<Output = T> + Copy,
{
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
