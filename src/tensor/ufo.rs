use crate::tensor::{CADJ, CAF, CAS, CF, CS, EUC, LOR, SPIN};

use super::{
    AbstractIndex, DenseTensor, HistoryStructure, IntoId,
    Representation::{self, Euclidean, Lorentz},
    SetTensorData, Shadowable, Slot, SparseTensor, TensorStructure, MAX_REP,
};
use num::{Complex, Float, One, Zero};

use symbolica::{
    representations::{Atom, Symbol},
    state::{State, Workspace},
};

pub const ID: Symbol = Symbol::init_var(MAX_REP + 1, 0);
pub const GAMMA: Symbol = Symbol::init_var(MAX_REP + 2, 0);
pub const GAMMA5: Symbol = Symbol::init_var(MAX_REP + 3, 0);
pub const PROJM: Symbol = Symbol::init_var(MAX_REP + 4, 0);
pub const PROJP: Symbol = Symbol::init_var(MAX_REP + 5, 0);
pub const SIGMA: Symbol = Symbol::init_var(MAX_REP + 6, 0);

pub fn init_state() {
    let mut state = State::get_global_state().write().unwrap();
    assert!(EUC == state.get_or_insert_fn("euc", None).unwrap());
    assert!(LOR == state.get_or_insert_fn("lor", None).unwrap());
    assert!(SPIN == state.get_or_insert_fn("spin", None).unwrap());
    assert!(CADJ == state.get_or_insert_fn("CAdj", None).unwrap());
    assert!(CF == state.get_or_insert_fn("CF", None).unwrap());
    assert!(CAF == state.get_or_insert_fn("CAF", None).unwrap());
    assert!(CS == state.get_or_insert_fn("CS", None).unwrap());
    assert!(CAS == state.get_or_insert_fn("CAS", None).unwrap());

    assert!(ID == state.get_or_insert_fn("id", None).unwrap());
    assert!(GAMMA == state.get_or_insert_fn("γ", None).unwrap());
    assert!(GAMMA5 == state.get_or_insert_fn("γ5", None).unwrap());
    assert!(PROJM == state.get_or_insert_fn("ProjM", None).unwrap());
    assert!(PROJP == state.get_or_insert_fn("ProjP", None).unwrap());
    assert!(SIGMA == state.get_or_insert_fn("σ", None).unwrap());
}

#[allow(dead_code)]
#[must_use]
pub fn identity<T>(
    indices: (AbstractIndex, AbstractIndex),
    signature: Representation,
) -> SparseTensor<Complex<T>>
where
    T: One + Zero,
{
    //TODO: make it just swap indices
    let structure = [(indices.0, signature), (indices.1, signature)]
        .into_iter()
        .map(Slot::from)
        .collect();
    let mut identity = SparseTensor::empty(structure);
    for i in 0..signature.into() {
        identity
            .set(&[i, i], Complex::<T>::new(T::one(), T::zero()))
            .unwrap_or_else(|_| unreachable!());
    }
    identity
}

/// Create a rank 2 identity tensor
///
/// # Arguments
///
/// * `structure` - The structure of the tensor
///
/// # Panics
///
/// * If the structure is not rank 2
/// * If the structure has different indices

pub fn identity_data<T, N>(structure: N) -> SparseTensor<T, N>
where
    T: One,
    N: TensorStructure,
{
    assert!(structure.order() == 2, "Identity tensor must be rank 2");

    assert!(
        structure.reps()[0] == structure.reps()[1],
        "Identity tensor must have equal indices"
    );

    let mut identity = SparseTensor::empty(structure);

    for i in 0..identity.shape()[0].into() {
        identity
            .set(&[i, i], T::one())
            .unwrap_or_else(|_| unreachable!());
    }
    identity
}

#[allow(dead_code)]
#[must_use]
pub fn lorentz_identity<T>(indices: (AbstractIndex, AbstractIndex)) -> SparseTensor<Complex<T>>
where
    T: One + Zero,
{
    // IdentityL(1,2) (Lorentz) Kronecker delta δ^μ1_μ1
    let signature = Lorentz(4.into());
    identity(indices, signature)
}

pub fn mink_four_vector<T>(index: AbstractIndex, p: &[T; 4]) -> DenseTensor<T>
where
    T: Clone,
{
    DenseTensor::from_data(p, vec![Slot::from((index, Lorentz(4.into())))])
        .unwrap_or_else(|_| unreachable!())
}

pub fn mink_four_vector_sym<T>(
    index: AbstractIndex,
    p: &[T; 4],
    state: &mut symbolica::state::State,
) -> DenseTensor<T, HistoryStructure<Symbol>>
where
    T: Clone,
{
    DenseTensor::from_data(
        p,
        HistoryStructure::new(
            &[(index, Lorentz(4.into()))],
            state
                .get_or_insert_fn("p", None)
                .unwrap_or_else(|_| unreachable!()),
        ),
    )
    .unwrap_or_else(|_| unreachable!())
}

pub fn euclidean_four_vector<T>(index: AbstractIndex, p: &[T; 4]) -> DenseTensor<T>
where
    T: Clone,
{
    DenseTensor::from_data(p, vec![Slot::from((index, Euclidean(4.into())))])
        .unwrap_or_else(|_| unreachable!())
}

pub fn euclidean_four_vector_sym<T>(
    index: AbstractIndex,
    p: &[T; 4],
    state: &mut symbolica::state::State,
) -> DenseTensor<T, HistoryStructure<Symbol>>
where
    T: Clone,
{
    DenseTensor::from_data(
        p,
        HistoryStructure::new(
            &[(index, Euclidean(4.into()))],
            state
                .get_or_insert_fn("p", None)
                .unwrap_or_else(|_| unreachable!()),
        ),
    )
    .unwrap_or_else(|_| unreachable!())
}

pub fn param_mink_four_vector<N>(
    index: AbstractIndex,
    name: N,
    _state: &mut symbolica::state::State,
    _ws: &Workspace,
) -> DenseTensor<Atom, HistoryStructure<N>>
where
    N: Clone + IntoId,
{
    HistoryStructure::new(&[(index, Lorentz(4.into()))], name)
        .shadow()
        .unwrap_or_else(|| unreachable!())
}

pub fn param_euclidean_four_vector<N>(
    index: AbstractIndex,
    name: N,
    _state: &mut symbolica::state::State,
    _ws: &Workspace,
) -> DenseTensor<Atom, HistoryStructure<N>>
where
    N: Clone + IntoId,
{
    HistoryStructure::new(&[(index, Euclidean(4.into()))], name)
        .shadow()
        .unwrap_or_else(|| unreachable!())
}

#[allow(dead_code)]
#[must_use]
pub fn euclidean_identity<T>(indices: (AbstractIndex, AbstractIndex)) -> SparseTensor<Complex<T>>
where
    T: One + Zero,
{
    // Identity(1,2) (Spinorial) Kronecker delta δ_s1_s2
    let signature = Euclidean(4.into());
    identity(indices, signature)
}

#[allow(dead_code)]
pub fn gamma<T>(
    minkindex: AbstractIndex,
    indices: (AbstractIndex, AbstractIndex),
) -> SparseTensor<Complex<T>>
where
    T: One + Zero + Copy + std::ops::Neg<Output = T>,
{
    // Gamma(1,2,3) Dirac matrix (γ^μ1)_s2_s3
    let structure = [
        (indices.0, Euclidean(4.into())),
        (indices.1, Euclidean(4.into())),
        (minkindex, Lorentz(4.into())),
    ]
    .into_iter()
    .map(Slot::from)
    .collect();

    gamma_data(structure)
}

pub fn gammasym<T>(
    minkindex: AbstractIndex,
    indices: (AbstractIndex, AbstractIndex),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, HistoryStructure<Symbol>>
where
    T: One + Zero + Copy + std::ops::Neg<Output = T>,
{
    let structure = HistoryStructure::new(
        &[
            (indices.0, Euclidean(4.into())),
            (indices.1, Euclidean(4.into())),
            (minkindex, Lorentz(4.into())),
        ],
        state
            .get_or_insert_fn("γ", None)
            .unwrap_or_else(|_| unreachable!()),
    );

    gamma_data(structure)
}

#[allow(clippy::similar_names)]
pub fn gamma_data<T, N>(structure: N) -> SparseTensor<Complex<T>, N>
where
    T: One + Zero + Copy + std::ops::Neg<Output = T>,
    N: TensorStructure,
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

pub fn gamma5<T>(indices: (AbstractIndex, AbstractIndex)) -> SparseTensor<Complex<T>>
where
    T: One + Zero + Copy,
{
    let structure = [
        (indices.0, Euclidean(4.into())),
        (indices.1, Euclidean(4.into())),
    ]
    .into_iter()
    .map(Slot::from)
    .collect();

    gamma5_data(structure)
}

pub fn gamma5sym<T>(
    indices: (AbstractIndex, AbstractIndex),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, HistoryStructure<Symbol>>
where
    T: One + Zero + Copy,
{
    let structure = HistoryStructure::new(
        &[
            (indices.0, Euclidean(4.into())),
            (indices.1, Euclidean(4.into())),
        ],
        state
            .get_or_insert_fn("γ5", None)
            .unwrap_or_else(|_| unreachable!()),
    );

    gamma5_data(structure)
}

pub fn gamma5_data<T, N>(structure: N) -> SparseTensor<Complex<T>, N>
where
    T: One + Zero + Copy,
    N: TensorStructure,
{
    let c1 = Complex::<T>::new(T::one(), T::zero());

    let mut gamma5 = SparseTensor::empty(structure);

    gamma5.set(&[0, 2], c1).unwrap();
    gamma5.set(&[1, 3], c1).unwrap();
    gamma5.set(&[2, 0], c1).unwrap();
    gamma5.set(&[3, 1], c1).unwrap();

    gamma5
}

pub fn proj_m<T>(indices: (AbstractIndex, AbstractIndex)) -> SparseTensor<Complex<T>>
where
    T: Float,
{
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
    let structure = [
        (indices.0, Euclidean(4.into())),
        (indices.1, Euclidean(4.into())),
    ]
    .into_iter()
    .map(Slot::from)
    .collect();

    proj_m_data(structure)
}

pub fn proj_msym<T>(
    indices: (AbstractIndex, AbstractIndex),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, HistoryStructure<Symbol>>
where
    T: Float,
{
    let structure = HistoryStructure::new(
        &[
            (indices.0, Euclidean(4.into())),
            (indices.1, Euclidean(4.into())),
        ],
        state
            .get_or_insert_fn("ProjM", None)
            .unwrap_or_else(|_| unreachable!()),
    );

    proj_m_data(structure)
}

#[allow(clippy::similar_names)]
pub fn proj_m_data<T, N>(structure: N) -> SparseTensor<Complex<T>, N>
where
    T: Float,
    N: TensorStructure,
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

pub fn proj_p<T>(indices: (AbstractIndex, AbstractIndex)) -> SparseTensor<Complex<T>>
where
    T: Float,
{
    // ProjP(1,2) Right chirality projector (( 1+γ5)/ 2 )_s1_s2
    let structure = [
        (indices.0, Euclidean(4.into())),
        (indices.1, Euclidean(4.into())),
    ]
    .into_iter()
    .map(Slot::from)
    .collect();

    proj_p_data(structure)
}

pub fn proj_psym<T>(
    indices: (AbstractIndex, AbstractIndex),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, HistoryStructure<Symbol>>
where
    T: Float,
{
    let structure = HistoryStructure::new(
        &[
            (indices.0, Euclidean(4.into())),
            (indices.1, Euclidean(4.into())),
        ],
        state
            .get_or_insert_fn("ProjP", None)
            .unwrap_or_else(|_| unreachable!()),
    );

    proj_p_data(structure)
}

pub fn proj_p_data<T, N>(structure: N) -> SparseTensor<Complex<T>, N>
where
    T: Float,
    N: TensorStructure,
{
    // ProjP(1,2) Right chirality projector (( 1+γ5)/ 2 )_s1_s2
    let chalf = Complex::<T>::new(T::from(0.5).unwrap_or_else(|| unreachable!()), T::zero());

    let mut proj_p = SparseTensor::empty(structure);

    proj_p
        .set(&[0, 0], chalf)
        .unwrap_or_else(|_| unreachable!());
    proj_p
        .set(&[1, 1], chalf)
        .unwrap_or_else(|_| unreachable!());
    proj_p
        .set(&[2, 2], chalf)
        .unwrap_or_else(|_| unreachable!());
    proj_p
        .set(&[3, 3], chalf)
        .unwrap_or_else(|_| unreachable!());

    proj_p
        .set(&[0, 2], chalf)
        .unwrap_or_else(|_| unreachable!());
    proj_p
        .set(&[1, 3], chalf)
        .unwrap_or_else(|_| unreachable!());
    proj_p
        .set(&[2, 0], chalf)
        .unwrap_or_else(|_| unreachable!());
    proj_p
        .set(&[3, 1], chalf)
        .unwrap_or_else(|_| unreachable!());

    proj_p
}

pub fn sigma<T>(
    indices: (AbstractIndex, AbstractIndex),
    minkdices: (AbstractIndex, AbstractIndex),
) -> SparseTensor<Complex<T>>
where
    T: One + Zero + std::ops::Neg<Output = T> + Copy,
{
    let structure = [
        (indices.0, Euclidean(4.into())),
        (indices.1, Euclidean(4.into())),
        (minkdices.0, Lorentz(4.into())),
        (minkdices.1, Lorentz(4.into())),
    ]
    .into_iter()
    .map(Slot::from)
    .collect();

    sigma_data(structure)
}

pub fn sigmasym<T>(
    indices: (AbstractIndex, AbstractIndex),
    minkdices: (AbstractIndex, AbstractIndex),
    state: &mut symbolica::state::State,
) -> SparseTensor<Complex<T>, HistoryStructure<Symbol>>
where
    T: One + Zero + std::ops::Neg<Output = T> + Copy,
{
    let structure = HistoryStructure::new(
        &[
            (indices.0, Euclidean(4.into())),
            (indices.1, Euclidean(4.into())),
            (minkdices.0, Lorentz(4.into())),
            (minkdices.1, Lorentz(4.into())),
        ],
        state
            .get_or_insert_fn("σ", None)
            .unwrap_or_else(|_| unreachable!()),
    );

    sigma_data(structure)
}

#[allow(clippy::similar_names)]
pub fn sigma_data<T, N>(structure: N) -> SparseTensor<Complex<T>, N>
where
    T: One + Zero + std::ops::Neg<Output = T> + Copy,
    N: TensorStructure,
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
