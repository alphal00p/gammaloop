use std::{ops::Neg, sync::LazyLock};

use idenso::{
    IndexTooling,
    color::CS,
    dirac::{AGS, spinor_matrix_structure},
    representations::initialize,
};

use spenso::{
    algebra::complex::Complex,
    network::{
        Network,
        library::{
            LibraryTensor, TensorLibraryData,
            function_lib::{INBUILTS, PanicMissingConcrete, SymbolLib},
            symbolic::{ExplicitKey, TensorLibrary},
        },
        parsing::ShadowedStructure,
        store::NetworkStore,
    },
    structure::{
        Canonicalized, TensorDataLayout, TensorStructure, abstract_index::AbstractIndex,
        slot::AbsInd,
    },
    tensors::{
        complex::RealOrComplexTensor,
        data::{SetTensorData, SparseTensor, StorageTensor},
        parametric::{MixedTensor, ParamTensor},
    },
};
use symbolica::{
    atom::{Atom, Symbol},
    parse_lit,
};

/// Nonzero Weyl-basis entries of the UFO charge-conjugation matrix C = -i γ² γ⁰.
pub const CHARGE_CONJUGATION_WEYL_COMPONENTS: [([usize; 2], i8); 4] =
    [([0, 1], -1), ([1, 0], 1), ([2, 3], 1), ([3, 2], -1)];

struct LogicalSparseInput<T, N> {
    layout: TensorDataLayout,
    storage: SparseTensor<T, N>,
}

impl<T, N: TensorStructure> LogicalSparseInput<T, N> {
    fn set(&mut self, indices: &[usize], value: T) -> Result<(), &'static str> {
        let storage_indices = self
            .layout
            .logical_expanded_to_storage_expanded(indices)
            .map_err(|_| "hard-coded logical coordinate must be in bounds")?;
        self.storage
            .set(&storage_indices, value)
            .map_err(|_| "translated storage coordinate must be in bounds")
    }
}

fn sparse_from_logical<T, N>(
    structure: Canonicalized<N>,
    zero: T,
    build: impl FnOnce(&mut LogicalSparseInput<T, N>),
) -> Canonicalized<SparseTensor<T, N>>
where
    N: TensorStructure,
{
    let layout = TensorDataLayout::from_canonicalized(&structure)
        .expect("hard-coded tensor data requires concrete dimensions");
    structure.map_canonical(|structure| {
        let mut input = LogicalSparseInput {
            layout,
            storage: SparseTensor::empty(structure, zero),
        };
        build(&mut input);
        input.storage
    })
}

#[allow(clippy::similar_names)]
pub fn gamma_data_dirac<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Clone + Neg<Output = T>,
    N: TensorStructure,
{
    let c1 = Complex::<T>::new(one.clone(), zero.clone());
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
    let ci = Complex::<T>::new(zero.clone(), one.clone());
    let cni = Complex::<T>::new(zero.clone(), -one.clone());
    sparse_from_logical(structure, z, |gamma| {
        // Coordinates follow the public logical order: bispinor row, bispinor
        // column, Minkowski. The layout maps them into canonical storage.

        // dirac gamma matrices

        gamma.set(&[0, 0, 0], c1.clone()).unwrap();
        gamma.set(&[1, 1, 0], c1.clone()).unwrap();
        gamma.set(&[2, 2, 0], cn1.clone()).unwrap();
        gamma.set(&[3, 3, 0], cn1.clone()).unwrap();

        gamma.set(&[0, 3, 1], c1.clone()).unwrap();
        gamma.set(&[1, 2, 1], c1.clone()).unwrap();
        gamma.set(&[2, 1, 1], cn1.clone()).unwrap();
        gamma.set(&[3, 0, 1], cn1.clone()).unwrap();

        gamma.set(&[0, 3, 2], cni.clone()).unwrap();
        gamma.set(&[1, 2, 2], ci.clone()).unwrap();
        gamma.set(&[2, 1, 2], ci.clone()).unwrap();
        gamma.set(&[3, 0, 2], cni.clone()).unwrap();

        gamma.set(&[0, 2, 3], c1.clone()).unwrap();
        gamma.set(&[1, 3, 3], cn1.clone()).unwrap();
        gamma.set(&[2, 0, 3], cn1.clone()).unwrap();
        gamma.set(&[3, 1, 3], c1.clone()).unwrap();

        // gamma.to_dense()
    })
}

#[allow(clippy::similar_names)]
pub fn gamma_data_weyl<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Neg<Output = T> + Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let c1 = Complex::<T>::new(one.clone(), zero.clone());
    let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
    let ci = Complex::<T>::new(zero.clone(), one.clone());
    let cni = Complex::<T>::new(zero.clone(), -one.clone());
    sparse_from_logical(structure, z, |gamma| {
        // Coordinates follow the public logical order: bispinor row, bispinor
        // column, Minkowski. The layout maps them into canonical storage.

        // dirac gamma matrices

        gamma.set(&[0, 2, 0], c1.clone()).unwrap();
        gamma.set(&[1, 3, 0], c1.clone()).unwrap();
        gamma.set(&[2, 0, 0], c1.clone()).unwrap();
        gamma.set(&[3, 1, 0], c1.clone()).unwrap();

        gamma.set(&[0, 3, 1], c1.clone()).unwrap();
        gamma.set(&[1, 2, 1], c1.clone()).unwrap();
        gamma.set(&[2, 1, 1], cn1.clone()).unwrap();
        gamma.set(&[3, 0, 1], cn1.clone()).unwrap();

        gamma.set(&[0, 3, 2], cni.clone()).unwrap();
        gamma.set(&[1, 2, 2], ci.clone()).unwrap();
        gamma.set(&[2, 1, 2], ci.clone()).unwrap();
        gamma.set(&[3, 0, 2], cni.clone()).unwrap();

        gamma.set(&[0, 2, 3], c1.clone()).unwrap();
        gamma.set(&[1, 3, 3], cn1.clone()).unwrap();
        gamma.set(&[2, 0, 3], cn1.clone()).unwrap();
        gamma.set(&[3, 1, 3], c1.clone()).unwrap();

        // gamma.to_dense()
    })
}

#[allow(clippy::similar_names)]
pub fn gamma_transpose_weyl<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Neg<Output = T> + Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let c1 = Complex::<T>::new(one.clone(), zero.clone());
    let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
    let ci = Complex::<T>::new(zero.clone(), one.clone());
    let cni = Complex::<T>::new(zero.clone(), -one.clone());
    sparse_from_logical(structure, z, |gamma| {
        // Coordinates follow the public logical order: bispinor row, bispinor
        // column, Minkowski. The layout maps them into canonical storage.

        // dirac gamma matrices

        gamma.set(&[2, 0, 0], c1.clone()).unwrap();
        gamma.set(&[3, 1, 0], c1.clone()).unwrap();
        gamma.set(&[0, 2, 0], c1.clone()).unwrap();
        gamma.set(&[1, 3, 0], c1.clone()).unwrap();

        gamma.set(&[3, 0, 1], c1.clone()).unwrap();
        gamma.set(&[2, 1, 1], c1.clone()).unwrap();
        gamma.set(&[1, 2, 1], cn1.clone()).unwrap();
        gamma.set(&[0, 3, 1], cn1.clone()).unwrap();

        gamma.set(&[3, 0, 2], cni.clone()).unwrap();
        gamma.set(&[2, 1, 2], ci.clone()).unwrap();
        gamma.set(&[1, 2, 2], ci.clone()).unwrap();
        gamma.set(&[0, 3, 2], cni.clone()).unwrap();

        gamma.set(&[2, 0, 3], c1.clone()).unwrap();
        gamma.set(&[3, 1, 3], cn1.clone()).unwrap();
        gamma.set(&[0, 2, 3], cn1.clone()).unwrap();
        gamma.set(&[1, 3, 3], c1.clone()).unwrap();

        // gamma.to_dense()
    })
}

#[allow(clippy::similar_names)]
pub fn gamma_conj_data_weyl<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Neg<Output = T> + Clone,
    N: TensorStructure,
{
    gamma_data_weyl(structure, one, zero).map_canonical(|tensor| {
        tensor.map_data(|a| {
            let Complex { re, im } = a;
            Complex { re, im: -im }
        })
    })
}

#[allow(clippy::similar_names)]
pub fn gamma_adj_data_weyl<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Neg<Output = T> + Clone,
    N: TensorStructure,
{
    gamma_transpose_weyl(structure, one, zero).map_canonical(|tensor| {
        tensor.map_data(|a| {
            let Complex { re, im } = a;
            Complex { re, im: -im }
        })
    })
}

pub fn gamma0_weyl<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Clone,
    N: TensorStructure,
{
    let c1 = Complex::<T>::new(one, zero.clone());
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    sparse_from_logical(structure, z, |gamma0| {
        // ! No check on actual structure, should expext bis,bis,lor

        // dirac gamma0 matrices

        gamma0.set(&[0, 2], c1.clone()).unwrap();
        gamma0.set(&[1, 3], c1.clone()).unwrap();
        gamma0.set(&[2, 0], c1.clone()).unwrap();
        gamma0.set(&[3, 1], c1.clone()).unwrap();
    })
}

pub fn gamma5_dirac_data<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Clone,
    N: TensorStructure,
{
    let c1 = Complex::<T>::new(one, zero.clone());

    let z = Complex::<T>::new(zero.clone(), zero.clone());
    sparse_from_logical(structure, z, |gamma5| {
        gamma5.set(&[0, 2], c1.clone()).unwrap();
        gamma5.set(&[1, 3], c1.clone()).unwrap();
        gamma5.set(&[2, 0], c1.clone()).unwrap();
        gamma5.set(&[3, 1], c1.clone()).unwrap();
    })
}

pub fn gamma5_weyl_data<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Clone + Neg<Output = T>,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let c1 = Complex::<T>::new(one, zero);

    sparse_from_logical(structure, z, |gamma5| {
        gamma5.set(&[0, 0], -c1.clone()).unwrap();
        gamma5.set(&[1, 1], -c1.clone()).unwrap();
        gamma5.set(&[2, 2], c1.clone()).unwrap();
        gamma5.set(&[3, 3], c1.clone()).unwrap();
    })
}

#[allow(clippy::similar_names)]
pub fn proj_m_data_dirac<T, N>(
    structure: Canonicalized<N>,
    half: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Clone + Neg<Output = T>,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2

    let chalf = Complex::<T>::new(half.clone(), zero.clone());
    let cnhalf = Complex::<T>::new(-half, zero);

    sparse_from_logical(structure, z, |proj_m| {
        proj_m.set(&[0, 0], chalf.clone()).unwrap();
        proj_m.set(&[1, 1], chalf.clone()).unwrap();
        proj_m.set(&[2, 2], chalf.clone()).unwrap();
        proj_m.set(&[3, 3], chalf.clone()).unwrap();

        proj_m.set(&[0, 2], cnhalf.clone()).unwrap();
        proj_m.set(&[1, 3], cnhalf.clone()).unwrap();
        proj_m.set(&[2, 0], cnhalf.clone()).unwrap();
        proj_m.set(&[3, 1], cnhalf.clone()).unwrap();
    })
}

#[allow(clippy::similar_names)]
pub fn proj_m_data_weyl<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
    let c1 = Complex::<T>::new(one, zero);
    sparse_from_logical(structure, z, |proj_m| {
        proj_m.set(&[0, 0], c1.clone()).unwrap();
        proj_m.set(&[1, 1], c1.clone()).unwrap();
    })
}

pub fn proj_p_data_dirac<T, N>(
    structure: Canonicalized<N>,
    half: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    // ProjP(1,2) Right chirality projector (( 1+γ5)/ 2 )_s1_s2
    let chalf = Complex::<T>::new(half, zero);

    sparse_from_logical(structure, z, |proj_p| {
        proj_p
            .set(&[0, 0], chalf.clone())
            .unwrap_or_else(|_| unreachable!());
        proj_p
            .set(&[1, 1], chalf.clone())
            .unwrap_or_else(|_| unreachable!());
        proj_p
            .set(&[2, 2], chalf.clone())
            .unwrap_or_else(|_| unreachable!());
        proj_p
            .set(&[3, 3], chalf.clone())
            .unwrap_or_else(|_| unreachable!());

        proj_p
            .set(&[0, 2], chalf.clone())
            .unwrap_or_else(|_| unreachable!());
        proj_p
            .set(&[1, 3], chalf.clone())
            .unwrap_or_else(|_| unreachable!());
        proj_p
            .set(&[2, 0], chalf.clone())
            .unwrap_or_else(|_| unreachable!());
        proj_p
            .set(&[3, 1], chalf.clone())
            .unwrap_or_else(|_| unreachable!());
    })
}

#[allow(clippy::similar_names)]
pub fn proj_p_data_weyl<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
    let c1 = Complex::<T>::new(one, zero);
    sparse_from_logical(structure, z, |proj_p| {
        proj_p.set(&[2, 2], c1.clone()).unwrap();
        proj_p.set(&[3, 3], c1.clone()).unwrap();
    })
}

/// Fundamental SU(3) generators in the normalization `Tr(T^a T^b)=1/2 delta^{ab}`.
///
/// The index order follows `CS.t_strct`: adjoint, fundamental, anti-fundamental.
pub fn su3_generator_data<N>(
    structure: Canonicalized<N>,
) -> Canonicalized<SparseTensor<Complex<f64>, N>>
where
    N: TensorStructure,
{
    let z = Complex::new(0., 0.);
    let h = 0.5;
    let s = 1. / (2. * 3_f64.sqrt());
    sparse_from_logical(structure, z, |t| {
        t.set(&[0, 0, 1], Complex::new(h, 0.)).unwrap();
        t.set(&[0, 1, 0], Complex::new(h, 0.)).unwrap();

        t.set(&[1, 0, 1], Complex::new(0., -h)).unwrap();
        t.set(&[1, 1, 0], Complex::new(0., h)).unwrap();

        t.set(&[2, 0, 0], Complex::new(h, 0.)).unwrap();
        t.set(&[2, 1, 1], Complex::new(-h, 0.)).unwrap();

        t.set(&[3, 0, 2], Complex::new(h, 0.)).unwrap();
        t.set(&[3, 2, 0], Complex::new(h, 0.)).unwrap();

        t.set(&[4, 0, 2], Complex::new(0., -h)).unwrap();
        t.set(&[4, 2, 0], Complex::new(0., h)).unwrap();

        t.set(&[5, 1, 2], Complex::new(h, 0.)).unwrap();
        t.set(&[5, 2, 1], Complex::new(h, 0.)).unwrap();

        t.set(&[6, 1, 2], Complex::new(0., -h)).unwrap();
        t.set(&[6, 2, 1], Complex::new(0., h)).unwrap();

        t.set(&[7, 0, 0], Complex::new(s, 0.)).unwrap();
        t.set(&[7, 1, 1], Complex::new(s, 0.)).unwrap();
        t.set(&[7, 2, 2], Complex::new(-2. * s, 0.)).unwrap();
    })
}

/// SU(3) structure constants for `[T^a,T^b]=i f^{abc} T^c`.
pub fn su3_structure_f_data<N>(structure: Canonicalized<N>) -> Canonicalized<SparseTensor<f64, N>>
where
    N: TensorStructure,
{
    fn set_antisymmetric<N>(
        f: &mut LogicalSparseInput<f64, N>,
        a: usize,
        b: usize,
        c: usize,
        value: f64,
    ) where
        N: TensorStructure,
    {
        f.set(&[a, b, c], value).unwrap();
        f.set(&[b, c, a], value).unwrap();
        f.set(&[c, a, b], value).unwrap();
        f.set(&[b, a, c], -value).unwrap();
        f.set(&[a, c, b], -value).unwrap();
        f.set(&[c, b, a], -value).unwrap();
    }

    sparse_from_logical(structure, 0., |f| {
        set_antisymmetric(f, 0, 1, 2, 1.);
        set_antisymmetric(f, 0, 3, 6, 0.5);
        set_antisymmetric(f, 0, 4, 5, -0.5);
        set_antisymmetric(f, 1, 3, 5, 0.5);
        set_antisymmetric(f, 1, 4, 6, 0.5);
        set_antisymmetric(f, 2, 3, 4, 0.5);
        set_antisymmetric(f, 2, 5, 6, -0.5);
        set_antisymmetric(f, 3, 4, 7, 3_f64.sqrt() / 2.);
        set_antisymmetric(f, 5, 6, 7, 3_f64.sqrt() / 2.);
    })
}

/// Exact Atom-backed SU(3) generators, useful for symbolic tensor libraries.
pub fn su3_generator_data_atom<N>(
    structure: Canonicalized<N>,
) -> Canonicalized<SparseTensor<Atom, N>>
where
    N: TensorStructure,
{
    let half = Atom::num(1) / Atom::num(2);
    let sqrt3 = parse_lit!(sqrt(3));
    let t8 = sqrt3.clone() / Atom::num(6);
    sparse_from_logical(structure, Atom::Zero, |t| {
        t.set(&[0, 0, 1], half.clone()).unwrap();
        t.set(&[0, 1, 0], half.clone()).unwrap();

        t.set(&[1, 0, 1], -Atom::i() * half.clone()).unwrap();
        t.set(&[1, 1, 0], Atom::i() * half.clone()).unwrap();

        t.set(&[2, 0, 0], half.clone()).unwrap();
        t.set(&[2, 1, 1], -half.clone()).unwrap();

        t.set(&[3, 0, 2], half.clone()).unwrap();
        t.set(&[3, 2, 0], half.clone()).unwrap();

        t.set(&[4, 0, 2], -Atom::i() * half.clone()).unwrap();
        t.set(&[4, 2, 0], Atom::i() * half.clone()).unwrap();

        t.set(&[5, 1, 2], half.clone()).unwrap();
        t.set(&[5, 2, 1], half.clone()).unwrap();

        t.set(&[6, 1, 2], -Atom::i() * half).unwrap();
        t.set(&[6, 2, 1], Atom::i() / Atom::num(2)).unwrap();

        t.set(&[7, 0, 0], t8.clone()).unwrap();
        t.set(&[7, 1, 1], t8.clone()).unwrap();
        t.set(&[7, 2, 2], -sqrt3 / Atom::num(3)).unwrap();
    })
}

/// Exact Atom-backed SU(3) structure constants.
pub fn su3_structure_f_data_atom<N>(
    structure: Canonicalized<N>,
) -> Canonicalized<SparseTensor<Atom, N>>
where
    N: TensorStructure,
{
    fn set_antisymmetric<N>(
        f: &mut LogicalSparseInput<Atom, N>,
        a: usize,
        b: usize,
        c: usize,
        value: Atom,
    ) where
        N: TensorStructure,
    {
        f.set(&[a, b, c], value.clone()).unwrap();
        f.set(&[b, c, a], value.clone()).unwrap();
        f.set(&[c, a, b], value.clone()).unwrap();
        f.set(&[b, a, c], -value.clone()).unwrap();
        f.set(&[a, c, b], -value.clone()).unwrap();
        f.set(&[c, b, a], -value).unwrap();
    }

    let half = Atom::num(1) / Atom::num(2);
    let sqrt3_half = parse_lit!(sqrt(3)) / Atom::num(2);

    sparse_from_logical(structure, Atom::Zero, |f| {
        set_antisymmetric(f, 0, 1, 2, Atom::num(1));
        set_antisymmetric(f, 0, 3, 6, half.clone());
        set_antisymmetric(f, 0, 4, 5, -half.clone());
        set_antisymmetric(f, 1, 3, 5, half.clone());
        set_antisymmetric(f, 1, 4, 6, half.clone());
        set_antisymmetric(f, 2, 3, 4, half.clone());
        set_antisymmetric(f, 2, 5, 6, -half);
        set_antisymmetric(f, 3, 4, 7, sqrt3_half.clone());
        set_antisymmetric(f, 5, 6, 7, sqrt3_half);
    })
}

#[allow(clippy::similar_names)]
pub fn sigma_data<T, N>(
    structure: Canonicalized<N>,
    one: T,
    zero: T,
) -> Canonicalized<SparseTensor<Complex<T>, N>>
where
    T: Clone + Neg<Output = T>,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let c1 = Complex::<T>::new(one.clone(), zero.clone());
    let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
    let ci = Complex::<T>::new(zero.clone(), one.clone());
    let cni = Complex::<T>::new(zero.clone(), -one.clone());

    sparse_from_logical(structure, z, |sigma| {
        sigma.set(&[0, 2, 0, 1], c1.clone()).unwrap();
        sigma.set(&[0, 2, 3, 0], c1.clone()).unwrap();
        sigma.set(&[0, 3, 1, 2], c1.clone()).unwrap();
        sigma.set(&[1, 0, 2, 2], c1.clone()).unwrap();
        sigma.set(&[1, 1, 1, 2], c1.clone()).unwrap();
        sigma.set(&[1, 3, 0, 2], c1.clone()).unwrap();
        sigma.set(&[2, 2, 1, 0], c1.clone()).unwrap();
        sigma.set(&[2, 2, 2, 1], c1.clone()).unwrap();
        sigma.set(&[2, 3, 3, 2], c1.clone()).unwrap();
        sigma.set(&[3, 0, 0, 2], c1.clone()).unwrap();
        sigma.set(&[3, 3, 2, 2], c1.clone()).unwrap();
        sigma.set(&[3, 1, 3, 2], c1.clone()).unwrap();
        sigma.set(&[0, 1, 3, 0], ci.clone()).unwrap();
        sigma.set(&[0, 3, 1, 1], ci.clone()).unwrap();
        sigma.set(&[0, 3, 2, 0], ci.clone()).unwrap();
        sigma.set(&[1, 0, 3, 3], ci.clone()).unwrap();
        sigma.set(&[1, 1, 0, 3], ci.clone()).unwrap();
        sigma.set(&[1, 1, 2, 0], ci.clone()).unwrap();
        sigma.set(&[2, 1, 1, 0], ci.clone()).unwrap();
        sigma.set(&[2, 3, 0, 0], ci.clone()).unwrap();
        sigma.set(&[2, 3, 3, 1], ci.clone()).unwrap();
        sigma.set(&[3, 0, 1, 3], ci.clone()).unwrap();
        sigma.set(&[3, 1, 0, 0], ci.clone()).unwrap();
        sigma.set(&[3, 1, 2, 3], ci.clone()).unwrap();
        sigma.set(&[0, 0, 3, 2], cn1.clone()).unwrap();
        sigma.set(&[0, 1, 0, 2], cn1.clone()).unwrap();
        sigma.set(&[0, 2, 1, 3], cn1.clone()).unwrap();
        sigma.set(&[1, 2, 0, 3], cn1.clone()).unwrap();
        sigma.set(&[1, 2, 1, 1], cn1.clone()).unwrap();
        sigma.set(&[1, 2, 2, 0], cn1.clone()).unwrap();
        sigma.set(&[2, 0, 1, 2], cn1.clone()).unwrap();
        sigma.set(&[2, 1, 2, 2], cn1.clone()).unwrap();
        sigma.set(&[2, 2, 3, 3], cn1.clone()).unwrap();
        sigma.set(&[3, 2, 0, 0], cn1.clone()).unwrap();
        sigma.set(&[3, 2, 2, 3], cn1.clone()).unwrap();
        sigma.set(&[3, 2, 3, 1], cn1.clone()).unwrap();
        sigma.set(&[0, 0, 2, 3], cni.clone()).unwrap();
        sigma.set(&[0, 0, 3, 1], cni.clone()).unwrap();
        sigma.set(&[0, 1, 1, 3], cni.clone()).unwrap();
        sigma.set(&[1, 0, 2, 1], cni.clone()).unwrap();
        sigma.set(&[1, 3, 0, 1], cni.clone()).unwrap();
        sigma.set(&[1, 3, 3, 0], cni.clone()).unwrap();
        sigma.set(&[2, 0, 0, 3], cni.clone()).unwrap();
        sigma.set(&[2, 0, 1, 1], cni.clone()).unwrap();
        sigma.set(&[2, 1, 3, 3], cni.clone()).unwrap();
        sigma.set(&[3, 0, 0, 1], cni.clone()).unwrap();
        sigma.set(&[3, 3, 1, 0], cni.clone()).unwrap();
        sigma.set(&[3, 3, 2, 1], cni.clone()).unwrap();
    })
}

pub fn hep_lib<Aind: AbsInd, T: TensorLibraryData + Clone + Default>(
    one: T,
    zero: T,
) -> TensorLibrary<MixedTensor<T, ExplicitKey<Aind>>, Aind>
where
{
    let mut weyl = TensorLibrary::new();
    initialize();
    weyl.update_ids();

    let gamma_key = gamma_data_weyl(AGS.gamma_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(Into::into);
    // println!("layout{:?}", gamma_key.layout());
    weyl.insert_explicit(gamma_key);
    let gamma_conj_key =
        gamma_conj_data_weyl(AGS.gamma_conj_strct::<Aind>(4), one.clone(), zero.clone())
            .map_canonical(Into::into);
    // println!("layout{:?}", gamma_key.layout());
    weyl.insert_explicit(gamma_conj_key);
    let gamma_adj_key =
        gamma_adj_data_weyl(AGS.gamma_adj_strct::<Aind>(4), one.clone(), zero.clone())
            .map_canonical(Into::into);
    // println!("layout{:?}", gamma_key.layout());
    weyl.insert_explicit(gamma_adj_key);
    let gamma0_key = gamma0_weyl(AGS.gamma0_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(Into::into);
    // println!("layout{:?}", gamma_key.layout());
    weyl.insert_explicit(gamma0_key);

    let gamma5_key = gamma5_weyl_data(AGS.gamma5_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(Into::into);
    weyl.insert_explicit(gamma5_key);

    let projm_key = proj_m_data_weyl(AGS.projm_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(Into::into);
    weyl.insert_explicit(projm_key);

    let projp_key = proj_p_data_weyl(AGS.projp_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(Into::into);
    weyl.insert_explicit(projp_key);

    let charge_conjugation = sparse_from_logical(
        spinor_matrix_structure::<Aind>(AGS.charge_conjugation, 4),
        Complex::new(zero.clone(), zero.clone()),
        |tensor| {
            for (indices, sign) in CHARGE_CONJUGATION_WEYL_COMPONENTS {
                let value = if sign < 0 { -one.clone() } else { one.clone() };
                tensor
                    .set(&indices, Complex::new(value, zero.clone()))
                    .unwrap();
            }
        },
    )
    .map_canonical(Into::into);
    weyl.insert_explicit(charge_conjugation);

    weyl
}

pub fn insert_su3_color_tensors<Aind: AbsInd>(
    lib: &mut TensorLibrary<MixedTensor<f64, ExplicitKey<Aind>>, Aind>,
) {
    initialize();

    let t_key = su3_generator_data(CS.t_strct::<Aind>(3, 8)).map_canonical(Into::into);
    lib.insert_explicit(t_key);

    let f_key = su3_structure_f_data(CS.f_strct::<Aind>(8)).map_canonical(Into::into);
    lib.insert_explicit(f_key);
}

pub fn hep_lib_su3<Aind: AbsInd>() -> TensorLibrary<MixedTensor<f64, ExplicitKey<Aind>>, Aind> {
    let mut lib = hep_lib(1., 0.);
    insert_su3_color_tensors(&mut lib);
    lib
}

pub fn hep_lib_atom<Aind: AbsInd, T>() -> TensorLibrary<T, Aind>
where
    T: LibraryTensor<Structure = ExplicitKey<Aind>>
        + SetTensorData<SetData = <T as LibraryTensor>::Data>
        + Clone
        + From<ParamTensor<ExplicitKey<Aind>>>,
    <T as LibraryTensor>::Data: TensorLibraryData,
{
    let mut weyl = TensorLibrary::new();
    initialize();
    weyl.update_ids_from::<ParamTensor<ExplicitKey<Aind>>>();

    let one = Atom::one();
    let zero = Atom::Zero;

    let gamma_key = gamma_data_weyl(AGS.gamma_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(|tensor| {
            ParamTensor::param(tensor.map_data(|a| a.re + a.im * Atom::i()).into()).into()
        });
    // println!("layout{:?}", gamma_key.layout());
    weyl.insert_explicit(gamma_key);
    let gamma_conj_key =
        gamma_conj_data_weyl(AGS.gamma_conj_strct::<Aind>(4), one.clone(), zero.clone())
            .map_canonical(|tensor| {
                ParamTensor::param(tensor.map_data(|a| a.re + a.im * Atom::i()).into()).into()
            });
    // println!("layout{:?}", gamma_key.layout());
    weyl.insert_explicit(gamma_conj_key);
    let gamma_adj_key =
        gamma_adj_data_weyl(AGS.gamma_adj_strct::<Aind>(4), one.clone(), zero.clone())
            .map_canonical(|tensor| {
                ParamTensor::param(tensor.map_data(|a| a.re + a.im * Atom::i()).into()).into()
            });
    // println!("layout{:?}", gamma_key.layout());
    weyl.insert_explicit(gamma_adj_key);
    let gamma0_key = gamma0_weyl(AGS.gamma0_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(|tensor| {
            ParamTensor::param(tensor.map_data(|a| a.re + a.im * Atom::i()).into()).into()
        });
    // println!("layout{:?}", gamma_key.layout());
    weyl.insert_explicit(gamma0_key);

    let gamma5_key = gamma5_weyl_data(AGS.gamma5_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(|tensor| {
            ParamTensor::param(tensor.map_data(|a| a.re + a.im * Atom::i()).into()).into()
        });
    weyl.insert_explicit(gamma5_key);

    let projm_key = proj_m_data_weyl(AGS.projm_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(|tensor| {
            ParamTensor::param(tensor.map_data(|a| a.re + a.im * Atom::i()).into()).into()
        });
    weyl.insert_explicit(projm_key);

    let projp_key = proj_p_data_weyl(AGS.projp_strct::<Aind>(4), one.clone(), zero.clone())
        .map_canonical(|tensor| {
            ParamTensor::param(tensor.map_data(|a| a.re + a.im * Atom::i()).into()).into()
        });
    weyl.insert_explicit(projp_key);

    let charge_conjugation = sparse_from_logical(
        spinor_matrix_structure::<Aind>(AGS.charge_conjugation, 4),
        Atom::Zero,
        |tensor| {
            for (indices, sign) in CHARGE_CONJUGATION_WEYL_COMPONENTS {
                tensor.set(&indices, Atom::num(sign)).unwrap();
            }
        },
    )
    .map_canonical(|tensor| ParamTensor::param(tensor.into()).into());
    weyl.insert_explicit(charge_conjugation);

    let color_t_key = su3_generator_data_atom(CS.t_strct::<Aind>(3, 8))
        .map_canonical(|tensor| ParamTensor::param(tensor.into()).into());
    weyl.insert_explicit(color_t_key);

    let color_f_key = su3_structure_f_data_atom(CS.f_strct::<Aind>(8))
        .map_canonical(|tensor| ParamTensor::param(tensor.into()).into());
    weyl.insert_explicit(color_f_key);

    weyl
}

pub type HepTensor<Aind> = MixedTensor<f64, ShadowedStructure<Aind>>;

pub type HepNet<Aind> =
    Network<NetworkStore<HepTensor<Aind>, Atom>, ExplicitKey<Aind>, Symbol, Aind>;

pub static HEP_LIB: LazyLock<
    TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>,
> = LazyLock::new(hep_lib_atom);

pub static FUN_LIB: LazyLock<
    SymbolLib<RealOrComplexTensor<f64, ShadowedStructure<AbstractIndex>>, PanicMissingConcrete>,
> = LazyLock::new(|| {
    let mut lib = PanicMissingConcrete::new_lib();
    lib.insert(INBUILTS.conj, |a| match a {
        RealOrComplexTensor::Complex(c) => RealOrComplexTensor::Complex(c.map_data(|x| x.conj())),
        RealOrComplexTensor::Real(r) => RealOrComplexTensor::Real(r),
    });
    lib.insert_scalar_fallible(INBUILTS.conj, |scalar| Ok(scalar.spenso_conj()));
    lib
});

#[cfg(test)]
mod tests {

    use idenso::representations::{Bispinor, ColorAdjoint, ColorFundamental};
    use spenso::{
        network::{
            ExecutionResult, MinIntermediateCost, Network, Sequential, SingleSmallestDegree,
            SmallestDegree, SmallestDegreeIter, Steps,
            library::symbolic::ETS,
            parsing::{ParseSettings, ShadowedStructure, StrictTensorFilter},
            store::NetworkStore,
            tags::SPENSO_TAG,
        },
        structure::{
            HasStructure,
            abstract_index::AbstractIndex,
            representation::{Minkowski, RepName},
        },
        tensors::data::{DenseTensor, GetTensorData, SparseOrDense},
    };
    use symbolica::{
        atom::{Atom, Symbol},
        parse, parse_lit,
    };

    use super::*;

    fn exact_default_scalar(expression: Atom) -> Atom {
        let mut network = HepNet::<AbstractIndex>::try_from_view(
            expression.as_view(),
            &*HEP_LIB,
            &ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps),
        )
        .unwrap();
        network
            .execute::<Sequential, MinIntermediateCost, _, _, _>(&*HEP_LIB, &*FUN_LIB)
            .unwrap();
        match network.result_tensor(&*HEP_LIB).unwrap() {
            ExecutionResult::One => Atom::one(),
            ExecutionResult::Zero => Atom::zero(),
            ExecutionResult::Val(value) => value.into_owned().scalar().unwrap().into(),
        }
    }

    #[test]
    fn default_library_keeps_generic_metrics_exact() {
        // Dimensions deliberately differ from the explicit four-dimensional HEP
        // entries: exactness belongs to the generic metric factories.
        let mink = Minkowski {}.new_rep(7);
        let bis = Bispinor {}.new_rep(6);
        let coad = ColorAdjoint {}.new_rep(11);
        let cof = ColorFundamental {}.new_rep(5);
        for (representations, dimension, lorentzian) in [
            ([mink.to_lib(); 2], 7, true),
            ([bis.to_lib(); 2], 6, false),
            ([coad.to_lib(); 2], 11, false),
            ([cof.to_lib(), cof.dual().to_lib()], 5, false),
        ] {
            let key = ExplicitKey::from_iter(representations, ETS.metric, None);
            let metric = HEP_LIB
                .get_storage(key.canonical())
                .unwrap()
                .into_owned()
                .try_into_parametric()
                .expect("exact default metric components")
                .to_dense();
            for row in 0..dimension {
                for column in 0..dimension {
                    let expected = if row != column {
                        0
                    } else if lorentzian && row > 0 {
                        -1
                    } else {
                        1
                    };
                    assert_eq!(
                        metric.get_owned([row, column]).unwrap(),
                        Atom::num(expected)
                    );
                }
            }
        }
        let numeric = hep_lib_su3::<AbstractIndex>();
        let key = ExplicitKey::from_iter([mink.to_lib(); 2], ETS.metric, None);
        assert!(
            numeric
                .get_storage(key.canonical())
                .unwrap()
                .into_owned()
                .try_into_concrete()
                .is_ok()
        );
    }

    #[test]
    fn default_component_contractions_preserve_exact_symbolic_coefficients() {
        initialize();
        let color = parse!(
            "x/7*f(coad(8,a),coad(8,b),coad(8,c))^2",
            default_namespace = "spenso"
        );
        assert_eq!(
            exact_default_scalar(color),
            parse!("24*x/7", default_namespace = "spenso")
        );
        let trace = parse!(
            "x/7*spenso::gamma(bis(4,i),bis(4,j),mink(4,mu))*spenso::gamma(bis(4,j),bis(4,i),mink(4,mu))",
            default_namespace = "spenso"
        );
        assert_eq!(
            exact_default_scalar(trace),
            parse!("16*x/7", default_namespace = "spenso")
        );

        let mut library =
            hep_lib_atom::<AbstractIndex, MixedTensor<f64, ExplicitKey<AbstractIndex>>>();
        let vector = SPENSO_TAG.rank_one_tensor_symbol("exact_default_momentum");
        let key = ExplicitKey::from_iter([Minkowski {}.new_rep(4).to_lib()], vector, None);
        library.insert_explicit(key.map_canonical(|structure| {
            MixedTensor::Param(ParamTensor::param(
                DenseTensor::from_storage_data(
                    vec![
                        parse!("x/3", default_namespace = "spenso"),
                        parse!("1/5"),
                        parse!("-2/7"),
                        parse!("3/11"),
                    ],
                    structure,
                )
                .unwrap()
                .into(),
            ))
        }));
        let source = parse!(
            "g(mink(4,a),mink(4,b))*exact_default_momentum(mink(4,a))*exact_default_momentum(mink(4,b))",
            default_namespace = "spenso"
        );
        let mut network =
            HepNet::try_from_view(source.as_view(), &library, &ParseSettings::default()).unwrap();
        network
            .execute::<Sequential, MinIntermediateCost, _, _, _>(&library, &*FUN_LIB)
            .unwrap();
        let ExecutionResult::Val(value) = network.result_tensor(&library).unwrap() else {
            panic!("nonzero momentum norm")
        };
        let value: Atom = value.into_owned().scalar().unwrap().into();
        assert_eq!(
            value,
            parse!("x^2/9-1/25-4/49-9/121", default_namespace = "spenso")
        );
    }

    #[test]
    fn default_components_annihilate_an_odd_color_network_exactly() {
        initialize();
        // K3,3 is odd under exchanging its first two left vertices: the three
        // right f tensors each exchange two ports. This proves zero without
        // using color-reduction formulas or the component contraction itself.
        let source = parse!(
            "g(coad(8,a),coad(8,j))*f(coad(8,j),coad(8,b),coad(8,c))*f(coad(8,d),coad(8,e),coad(8,f))*f(coad(8,g),coad(8,h),coad(8,i))*f(coad(8,a),coad(8,d),coad(8,g))*f(coad(8,b),coad(8,e),coad(8,h))*f(coad(8,c),coad(8,f),coad(8,i))",
            default_namespace = "spenso"
        );
        assert!(!source.is_zero());
        assert_eq!(exact_default_scalar(source), Atom::zero());
    }

    #[test]
    fn su3_color_traces_are_independent_of_contraction_order() {
        use spenso::network::{MinIntermediateCost, MinResultRank, library::function_lib::Wrap};

        initialize();
        type Net = Network<
            NetworkStore<ParamTensor<ShadowedStructure<AbstractIndex>>, Atom>,
            ExplicitKey<AbstractIndex>,
            Symbol,
        >;
        let lib: TensorLibrary<ParamTensor<ExplicitKey<AbstractIndex>>, AbstractIndex> =
            hep_lib_atom();
        let functions: SymbolLib<ParamTensor<ShadowedStructure<AbstractIndex>>, Wrap> =
            Wrap::new_lib();
        for (expression, expected) in [
            // The two crossed cyclic orderings used to give -1/6 rather than -2/3.
            (
                "g(coad(8,a),coad(8,b))*g(coad(8,c),coad(8,d))*trace(cof(3),t(coad(8,a),in,out),t(coad(8,c),in,out),t(coad(8,b),in,out),t(coad(8,d),in,out))",
                parse!("-2/3"),
            ),
            (
                "1/8*(2+g(coad(8,a),coad(8,b))*g(coad(8,c),coad(8,d))*trace(cof(3),sym(t(coad(8,a),in,out),t(coad(8,b),in,out),t(coad(8,c),in,out),t(coad(8,d),in,out))))",
                parse!("2/3"),
            ),
            // sum_d T_d X T_d = Tr(X) I/2 - X/6 and
            // f_abc Tr(T_a T_b T_c) = 6i independently give -i.
            (
                "f(coad(8,a),coad(8,b),coad(8,c))*t(coad(8,a),cof(3,i0),dind(cof(3,i1)))*t(coad(8,d),cof(3,i1),dind(cof(3,i2)))*t(coad(8,b),cof(3,i2),dind(cof(3,i3)))*t(coad(8,c),cof(3,i3),dind(cof(3,i4)))*t(coad(8,d),cof(3,i4),dind(cof(3,i0)))",
                parse!("-1i"),
            ),
        ] {
            let expression = parse!(expression, default_namespace = "spenso");
            let original =
                Net::try_from_view(expression.as_view(), &lib, &ParseSettings::default()).unwrap();
            for prepare in [false, true] {
                let mut prepared = original.clone();
                if prepare {
                    prepared.graph.contract_ready_sum_boundaries();
                }
                macro_rules! check {
                    ($strategy:ty) => {{
                        let mut network = prepared.clone();
                        network
                            .execute::<Sequential, $strategy, _, _, _>(&lib, &functions)
                            .unwrap();
                        let ExecutionResult::Val(actual) = network.result_scalar().unwrap() else {
                            panic!("expected a scalar color trace");
                        };
                        assert_eq!(
                            actual.as_ref(),
                            &expected,
                            "strategy={}, prepare={prepare}",
                            stringify!($strategy)
                        );
                    }};
                }
                check!(SmallestDegree);
                check!(MinIntermediateCost);
                check!(MinResultRank);
            }
        }
    }

    #[test]
    fn scalar_conjugation_execution_preserves_exact_coefficients() {
        initialize();
        let mut network =
            HepNet::<AbstractIndex>::from_scalar(parse!("1/3 + 2i/7")).fun(INBUILTS.conj);
        network
            .execute::<Sequential, SmallestDegree, _, _, _>(&*HEP_LIB, &*FUN_LIB)
            .unwrap();
        let ExecutionResult::Val(conjugated) = network.result_scalar().unwrap() else {
            panic!("conjugation should produce a scalar value")
        };

        assert_eq!(conjugated.into_owned(), parse!("1/3 - 2i/7"));
    }

    #[test]
    fn dirac_gamma_data_uses_storage_order() {
        initialize();
        let gamma =
            gamma_data_dirac(AGS.gamma_strct::<AbstractIndex>(4), 1_i32, 0_i32).into_canonical();

        // The final coordinate is the Minkowski component. These are the
        // diagonal entries of gamma^0 in the Dirac basis.
        assert_eq!(*gamma.get_ref([0, 0, 0]).unwrap(), Complex::new(1, 0));
        assert_eq!(*gamma.get_ref([1, 1, 0]).unwrap(), Complex::new(1, 0));
        assert_eq!(*gamma.get_ref([2, 2, 0]).unwrap(), Complex::new(-1, 0));
        assert_eq!(*gamma.get_ref([3, 3, 0]).unwrap(), Complex::new(-1, 0));
    }

    #[test]
    fn charge_conjugation_weyl_components_and_clifford_identities() {
        initialize();
        let key = spinor_matrix_structure::<AbstractIndex>(AGS.charge_conjugation, 4);
        let concrete_library = hep_lib::<AbstractIndex, i32>(1, 0);
        let concrete = concrete_library
            .get_storage(key.canonical())
            .unwrap()
            .into_owned()
            .try_into_concrete()
            .unwrap();
        let RealOrComplexTensor::Complex(charge_conjugation) = concrete else {
            panic!("the Weyl library stores complex matrix components");
        };
        let charge_conjugation = charge_conjugation.to_dense();
        let atom_library =
            hep_lib_atom::<AbstractIndex, MixedTensor<i32, ExplicitKey<AbstractIndex>>>();
        let atom_matrix = atom_library
            .get_storage(key.canonical())
            .unwrap()
            .into_owned()
            .try_into_parametric()
            .unwrap()
            .to_dense();
        let gamma = gamma_data_weyl(AGS.gamma_strct::<AbstractIndex>(4), 1, 0)
            .into_canonical()
            .to_dense();

        for row in 0..4 {
            for column in 0..4 {
                let value = charge_conjugation.get_owned([row, column]).unwrap();
                // Derive the matrix independently from the shared Weyl gammas.
                let definition = (0..4).fold(Complex::new(0, 0), |sum, inner| {
                    sum + Complex::new(0, -1)
                        * gamma.get_owned([row, inner, 2]).unwrap()
                        * gamma.get_owned([inner, column, 0]).unwrap()
                });
                assert_eq!(value, definition);
                assert_eq!(value.conj(), value);
                assert_eq!(value, -charge_conjugation.get_owned([column, row]).unwrap());
                assert_eq!(
                    atom_matrix.get_owned([row, column]).unwrap(),
                    Atom::num(value.re)
                );

                let square = (0..4).fold(Complex::new(0, 0), |sum, inner| {
                    sum + charge_conjugation.get_owned([row, inner]).unwrap()
                        * charge_conjugation.get_owned([inner, column]).unwrap()
                });
                assert_eq!(square, Complex::new(-i32::from(row == column), 0));
                for mu in 0..4 {
                    let mut sandwich = Complex::new(0, 0);
                    let mut transposed_sandwich = Complex::new(0, 0);
                    for first in 0..4 {
                        for second in 0..4 {
                            let left = charge_conjugation.get_owned([row, first]).unwrap();
                            let right = charge_conjugation.get_owned([second, column]).unwrap();
                            sandwich +=
                                left * gamma.get_owned([first, second, mu]).unwrap() * right;
                            transposed_sandwich +=
                                left * gamma.get_owned([second, first, mu]).unwrap() * right;
                        }
                    }
                    assert_eq!(sandwich, gamma.get_owned([column, row, mu]).unwrap());
                    assert_eq!(
                        transposed_sandwich,
                        gamma.get_owned([row, column, mu]).unwrap()
                    );
                }
            }
        }
    }

    #[test]
    fn simple_scalar() {
        initialize();
        let gamma = AGS.gamma_strct(4);
        let _a = HEP_LIB.get(gamma.canonical()).unwrap();

        let expr = parse!("gamma(bis(4,l_5),bis(4,l_4),mink(4,l_4))*gamma(bis(4,l_6),bis(4,l_5),mink(4,l_4))*gamma(bis(4,l_4),bis(4,l_6),mink(4,l_5))*p(mink(4,l_5))
            ",default_namespace="spenso");
        // let expr = parse!(
        // "gamma(bis(4,l_4),bis(4,l_6),mink(4,l_5))*p(mink(4,l_5))
        // ",
        // "spenso"
        // );
        // println!("{}", expr);

        let mut net = Network::<
            NetworkStore<MixedTensor<f64, ShadowedStructure<AbstractIndex>>, Atom>,
            _,
            Symbol,
        >::try_from_view(
            expr.as_view(),
            &*HEP_LIB,
            &ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps),
        )
        .unwrap();

        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |a| Some(format!("{}", a.global_name.unwrap())),
                |a| a.structure().global_name.unwrap().to_string(),
                |a| a.to_string()
            )
        );

        net.execute::<Steps<1>, SmallestDegreeIter<1>, _, _, _>(&*HEP_LIB, &*FUN_LIB)
            .unwrap();
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |a| Some(format!("{}", a.global_name?)),
                |a| a
                    .structure()
                    .global_name
                    .map(|a| a.to_string())
                    .unwrap_or("".to_string()),
                |a| a.to_string()
            )
        );
        net.execute::<Steps<1>, SmallestDegreeIter<2>, _, _, _>(&*HEP_LIB, &*FUN_LIB)
            .unwrap();
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |a| Some(format!("{}", a.global_name?)),
                |a| a
                    .structure()
                    .global_name
                    .map(|a| a.to_string())
                    .unwrap_or("".to_string()),
                |a| a.to_string()
            )
        );

        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.to_string(),
                |a| a.to_string()
            )
        );
        // if let ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) =
        //     net.result().unwrap()
        // {
        //     // println!("YaY:{}", (&expr - &tensor.expression).expand());
        //     // assert_eq!(expr, tensor.expression);
        // } else {
        //     panic!("Not tensor")
        // }
    }

    #[test]
    // #[should_panic]
    fn parse_problem() {
        initialize();
        let gamma = AGS.gamma_strct(4);
        let _a = HEP_LIB.get(gamma.canonical()).unwrap();

        let expr = parse_lit!(
            (-1 * G
                ^ 3 * P(0, mink(4, 0))
                    * P(2, mink(4, 26))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 1))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    + -1 * G
                ^ 3 * P(0, mink(4, 26))
                    * P(1, mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 0))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    + -1 * G
                ^ 3 * P(0, mink(4, 26))
                    * P(1, mink(4, 5))
                    * g(mink(4, 0), mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 5))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    + -1 * G
                ^ 3 * P(0, mink(4, 5))
                    * P(2, mink(4, 26))
                    * g(mink(4, 0), mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 5))
                    + -1 * G
                ^ 3 * P(1, mink(4, 1))
                    * P(1, mink(4, 26))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 0))
                    + -1 * G
                ^ 3 * P(1, mink(4, 26))
                    * P(1, mink(4, 5))
                    * g(mink(4, 0), mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 5))
                    + -2 * G
                ^ 3 * P(0, mink(4, 1))
                    * P(0, mink(4, 26))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 0))
                    + -2 * G
                ^ 3 * P(0, mink(4, 1))
                    * P(1, mink(4, 26))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 0))
                    + -2 * G
                ^ 3 * P(0, mink(4, 5))
                    * Q(0, mink(4, 5))
                    * g(mink(4, 0), mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 2), mink(4, 4))
                    + -2 * G
                ^ 3 * P(1, mink(4, 0))
                    * P(1, mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 2), mink(4, 4))
                    + -2 * G
                ^ 3 * P(1, mink(4, 0))
                    * P(2, mink(4, 26))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 1))
                    + -2 * G
                ^ 3 * P(1, mink(4, 1))
                    * P(2, mink(4, 0))
                    * gamma(bis(4, 3), bis(4, 2), mink(4, 4))
                    + -2 * G
                ^ 3 * P(1, mink(4, 5))
                    * P(2, mink(4, 5))
                    * g(mink(4, 0), mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 2), mink(4, 4))
                    + -4 * G
                ^ 3 * P(0, mink(4, 1))
                    * P(2, mink(4, 0))
                    * gamma(bis(4, 3), bis(4, 2), mink(4, 4))
                    + 2 * G
                ^ 3 * P(0, mink(4, 0))
                    * P(0, mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 2), mink(4, 4))
                    + 2 * G
                ^ 3 * P(0, mink(4, 0))
                    * P(2, mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 2), mink(4, 4))
                    + 2 * G
                ^ 3 * P(0, mink(4, 1))
                    * P(2, mink(4, 26))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 0))
                    + 2 * G
                ^ 3 * P(0, mink(4, 26))
                    * P(1, mink(4, 0))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 1))
                    + 2 * G
                ^ 3 * P(0, mink(4, 5))
                    * P(2, mink(4, 5))
                    * g(mink(4, 0), mink(4, 1))
                    * gamma(bis(4, 3), bis(4, 2), mink(4, 4))
                    + 2 * G
                ^ 3 * P(1, mink(4, 0))
                    * P(1, mink(4, 26))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 1))
                    + 2 * G
                ^ 3 * P(1, mink(4, i))
                ^ 2 * g(mink(4, 0), mink(4, 1)) * gamma(bis(4, 3), bis(4, 2), mink(4, 4)) + G
                ^ 3 * P(1, mink(4, 1))
                    * P(2, mink(4, 26))
                    * gamma(bis(4, 3), bis(4, 7), mink(4, 4))
                    * gamma(bis(4, 6), bis(4, 2), mink(4, 26))
                    * gamma(bis(4, 7), bis(4, 6), mink(4, 0))),
            default_namespace = "spenso"
        );
        // println!("{}", expr);

        let mut net = Network::<
            NetworkStore<MixedTensor<f64, ShadowedStructure<AbstractIndex>>, Atom>,
            _,
            Symbol,
        >::try_from_view(
            expr.as_view(),
            &*HEP_LIB,
            &ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps),
        )
        .unwrap();

        net.merge_ops();
        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |a| Some(format!("{}", a.global_name.unwrap())),
                |a| a.structure().global_name.unwrap().to_string(),
                |a| a.to_string()
            )
        );

        // net.validate();
        net.execute::<Steps<1>, SingleSmallestDegree<true>, _, _, _>(&*HEP_LIB, &(*FUN_LIB))
            .unwrap();
        // net.validate();
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.validate();
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.validate();
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.validate();
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.validate();
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.validate();
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.validate();
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.validate();
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.validate();
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.validate();
        // net.execute::<Steps<1>, ContractScalars, _, _>(&*HEP_LIB);
        // net.execute::<Steps<1>, SmallestDegree, _, _>(&*HEP_LIB);
        // net.execute::<StepsDebug<1>, SingleSmallestDegree<true>, _, _>(&*HEP_LIB);
        // net.execute::<Steps<1>, ContractScalars, _, _>(&*HEP_LIB);

        //     .unwrap();
        // net.execute::<Steps<14>, SingleSmallestDegree<false>, _, _>(&*HEP_LIB)
        //     .unwrap();
        // net.execute::<Steps<1>, SingleSmallestDegree<true>, _, _>(&*HEP_LIB)
        //     .unwrap();
        // // net.execute::<Sequential, SmallestDegree, _, _>(&*HEP_LIB)
        //     .unwrap();
        // println!(
        //     "{}",
        //     net.dot_display_impl(|a| a.to_string(), |_| None, |a| a.to_string())
        // );

        println!(
            "{}",
            net.dot_display_impl(
                |a| a.to_string(),
                |_| None,
                |a| a.structure().to_string().replace('\n', "\\n"),
                |a| a.to_string()
            )
        );
        // if let ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) =
        //     net.result().unwrap()
        // {
        //     // println!("YaY:{}", (&expr - &tensor.expression).expand());
        //     // assert_eq!(expr, tensor.expression);
        // } else {
        //     panic!("Not tensor")
        // }
    }

    #[test]
    fn transpose_test() {}
}
