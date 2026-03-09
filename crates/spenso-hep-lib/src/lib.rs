use std::{ops::Neg, sync::LazyLock};

use idenso::{gamma::AGS, representations::initialize};

use spenso::{
    algebra::complex::Complex,
    network::{
        Network,
        library::{
            TensorLibraryData,
            function_lib::{INBUILTS, PanicMissingConcrete, SymbolLib},
            symbolic::{ExplicitKey, TensorLibrary},
        },
        parsing::ShadowedStructure,
        store::NetworkStore,
    },
    structure::{PermutedStructure, TensorStructure, abstract_index::AbstractIndex, slot::AbsInd},
    tensors::{
        complex::RealOrComplexTensor,
        data::{SetTensorData, SparseTensor, StorageTensor},
        parametric::{MixedTensor, ParamOrConcrete},
    },
};
use symbolica::atom::{Atom, Symbol};

#[allow(clippy::similar_names)]
pub fn gamma_data_dirac<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Clone + Neg<Output = T>,
    N: TensorStructure,
{
    let c1 = Complex::<T>::new(one.clone(), zero.clone());
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
    let ci = Complex::<T>::new(zero.clone(), one.clone());
    let cni = Complex::<T>::new(zero.clone(), -one.clone());
    let mut gamma = SparseTensor::empty(structure, z);
    // ! No check on actual structure, should expext mink,bis,bis

    // dirac gamma matrices

    gamma.set(&[0, 0, 0], c1.clone()).unwrap();
    gamma.set(&[0, 1, 1], c1.clone()).unwrap();
    gamma.set(&[0, 2, 2], cn1.clone()).unwrap();
    gamma.set(&[0, 3, 3], cn1.clone()).unwrap();

    gamma.set(&[1, 0, 3], c1.clone()).unwrap();
    gamma.set(&[1, 1, 2], c1.clone()).unwrap();
    gamma.set(&[1, 2, 1], cn1.clone()).unwrap();
    gamma.set(&[1, 3, 0], cn1.clone()).unwrap();

    gamma.set(&[2, 0, 3], cni.clone()).unwrap();
    gamma.set(&[2, 1, 2], ci.clone()).unwrap();
    gamma.set(&[2, 2, 1], ci.clone()).unwrap();
    gamma.set(&[2, 3, 0], cni.clone()).unwrap();

    gamma.set(&[3, 0, 2], c1.clone()).unwrap();
    gamma.set(&[3, 1, 3], cn1.clone()).unwrap();
    gamma.set(&[3, 2, 0], cn1.clone()).unwrap();
    gamma.set(&[3, 3, 1], c1.clone()).unwrap();

    gamma //.to_dense()
}

#[allow(clippy::similar_names)]
pub fn gamma_data_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Neg<Output = T> + Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let c1 = Complex::<T>::new(one.clone(), zero.clone());
    let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
    let ci = Complex::<T>::new(zero.clone(), one.clone());
    let cni = Complex::<T>::new(zero.clone(), -one.clone());
    let mut gamma = SparseTensor::empty(structure, z);
    // ! No check on actual structure, should expext mink,bis,bis

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

    gamma //.to_dense()
}

#[allow(clippy::similar_names)]
pub fn gamma_transpose_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Neg<Output = T> + Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let c1 = Complex::<T>::new(one.clone(), zero.clone());
    let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
    let ci = Complex::<T>::new(zero.clone(), one.clone());
    let cni = Complex::<T>::new(zero.clone(), -one.clone());
    let mut gamma = SparseTensor::empty(structure, z);
    // ! No check on actual structure, should expext mink,bis,bis

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

    gamma //.to_dense()
}

#[allow(clippy::similar_names)]
pub fn gamma_conj_data_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Neg<Output = T> + Clone,
    N: TensorStructure,
{
    gamma_data_weyl(structure, one, zero).map_data(|a| {
        let Complex { re, im } = a;
        Complex { re, im: -im }
    })
}

#[allow(clippy::similar_names)]
pub fn gamma_adj_data_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Neg<Output = T> + Clone,
    N: TensorStructure,
{
    gamma_transpose_weyl(structure, one, zero).map_data(|a| {
        let Complex { re, im } = a;
        Complex { re, im: -im }
    })
}

pub fn gamma0_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Clone,
    N: TensorStructure,
{
    let c1 = Complex::<T>::new(one, zero.clone());
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let mut gamma0 = SparseTensor::empty(structure, z);
    // ! No check on actual structure, should expext bis,bis,lor

    // dirac gamma0 matrices

    gamma0.set(&[0, 2], c1.clone()).unwrap();
    gamma0.set(&[1, 3], c1.clone()).unwrap();
    gamma0.set(&[2, 0], c1.clone()).unwrap();
    gamma0.set(&[3, 1], c1.clone()).unwrap();

    gamma0
}

pub fn gamma5_dirac_data<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Clone,
    N: TensorStructure,
{
    let c1 = Complex::<T>::new(one, zero.clone());

    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let mut gamma5 = SparseTensor::empty(structure, z);

    gamma5.set(&[0, 2], c1.clone()).unwrap();
    gamma5.set(&[1, 3], c1.clone()).unwrap();
    gamma5.set(&[2, 0], c1.clone()).unwrap();
    gamma5.set(&[3, 1], c1.clone()).unwrap();

    gamma5
}

pub fn gamma5_weyl_data<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Clone + Neg<Output = T>,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let c1 = Complex::<T>::new(one, zero);

    let mut gamma5 = SparseTensor::empty(structure, z);

    gamma5.set(&[0, 0], -c1.clone()).unwrap();
    gamma5.set(&[1, 1], -c1.clone()).unwrap();
    gamma5.set(&[2, 2], c1.clone()).unwrap();
    gamma5.set(&[3, 3], c1.clone()).unwrap();

    gamma5
}

#[allow(clippy::similar_names)]
pub fn proj_m_data_dirac<T, N>(structure: N, half: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Clone + Neg<Output = T>,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2

    let chalf = Complex::<T>::new(half.clone(), zero.clone());
    let cnhalf = Complex::<T>::new(-half, zero);

    let mut proj_m = SparseTensor::empty(structure, z);

    proj_m.set(&[0, 0], chalf.clone()).unwrap();
    proj_m.set(&[1, 1], chalf.clone()).unwrap();
    proj_m.set(&[2, 2], chalf.clone()).unwrap();
    proj_m.set(&[3, 3], chalf.clone()).unwrap();

    proj_m.set(&[0, 2], cnhalf.clone()).unwrap();
    proj_m.set(&[1, 3], cnhalf.clone()).unwrap();
    proj_m.set(&[2, 0], cnhalf.clone()).unwrap();
    proj_m.set(&[3, 1], cnhalf.clone()).unwrap();

    proj_m
}

#[allow(clippy::similar_names)]
pub fn proj_m_data_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
    let c1 = Complex::<T>::new(one, zero);
    let mut proj_m = SparseTensor::empty(structure, z);

    proj_m.set(&[0, 0], c1.clone()).unwrap();
    proj_m.set(&[1, 1], c1.clone()).unwrap();

    proj_m
}

pub fn proj_p_data_dirac<T, N>(structure: N, half: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    // ProjP(1,2) Right chirality projector (( 1+γ5)/ 2 )_s1_s2
    let chalf = Complex::<T>::new(half, zero);

    let mut proj_p = SparseTensor::empty(structure, z);

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

    proj_p
}

#[allow(clippy::similar_names)]
pub fn proj_p_data_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Clone,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
    let c1 = Complex::<T>::new(one, zero);
    let mut proj_p = SparseTensor::empty(structure, z);

    proj_p.set(&[2, 2], c1.clone()).unwrap();
    proj_p.set(&[3, 3], c1.clone()).unwrap();

    proj_p
}

#[allow(clippy::similar_names)]
pub fn sigma_data<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
where
    T: Clone + Neg<Output = T>,
    N: TensorStructure,
{
    let z = Complex::<T>::new(zero.clone(), zero.clone());
    let c1 = Complex::<T>::new(one.clone(), zero.clone());
    let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
    let ci = Complex::<T>::new(zero.clone(), one.clone());
    let cni = Complex::<T>::new(zero.clone(), -one.clone());

    let mut sigma = SparseTensor::empty(structure, z);
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

    sigma
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

    let gamma_key = PermutedStructure::identity(
        gamma_data_weyl(AGS.gamma_strct::<Aind>(4), one.clone(), zero.clone()).into(),
    );
    // println!("permutation{}", gamma_key.rep_permutation);
    weyl.insert_explicit(gamma_key);
    let gamma_conj_key = PermutedStructure::identity(
        gamma_conj_data_weyl(AGS.gamma_conj_strct::<Aind>(4), one.clone(), zero.clone()).into(),
    );
    // println!("permutation{}", gamma_key.rep_permutation);
    weyl.insert_explicit(gamma_conj_key);
    let gamma_adj_key = PermutedStructure::identity(
        gamma_adj_data_weyl(AGS.gamma_adj_strct::<Aind>(4), one.clone(), zero.clone()).into(),
    );
    // println!("permutation{}", gamma_key.rep_permutation);
    weyl.insert_explicit(gamma_adj_key);
    let gamma0_key = PermutedStructure::identity(
        gamma0_weyl(AGS.gamma0_strct::<Aind>(4), one.clone(), zero.clone()).into(),
    );
    // println!("permutation{}", gamma_key.rep_permutation);
    weyl.insert_explicit(gamma0_key);

    let gamma5_key = PermutedStructure::identity(
        gamma5_weyl_data(AGS.gamma5_strct::<Aind>(4), one.clone(), zero.clone()).into(),
    );
    weyl.insert_explicit(gamma5_key);

    let projm_key = PermutedStructure::identity(
        proj_m_data_weyl(AGS.projm_strct::<Aind>(4), one.clone(), zero.clone()).into(),
    );
    weyl.insert_explicit(projm_key);

    let projp_key = PermutedStructure::identity(
        proj_p_data_weyl(AGS.projp_strct::<Aind>(4), one.clone(), zero.clone()).into(),
    );
    weyl.insert_explicit(projp_key);

    weyl
}

pub fn hep_lib_atom<Aind: AbsInd, T: TensorLibraryData + Clone + Default>()
-> TensorLibrary<MixedTensor<T, ExplicitKey<Aind>>, Aind>
where
{
    let mut weyl = TensorLibrary::new();
    initialize();
    weyl.update_ids();

    let one = Atom::one();
    let zero = Atom::Zero;

    let gamma_key = PermutedStructure::identity(ParamOrConcrete::param(
        gamma_data_weyl(AGS.gamma_strct::<Aind>(4), one.clone(), zero.clone())
            .map_data(|a| a.re + a.im * Atom::i())
            .into(),
    ));
    // println!("permutation{}", gamma_key.rep_permutation);
    weyl.insert_explicit(gamma_key);
    let gamma_conj_key = PermutedStructure::identity(ParamOrConcrete::param(
        gamma_conj_data_weyl(AGS.gamma_conj_strct::<Aind>(4), one.clone(), zero.clone())
            .map_data(|a| a.re + a.im * Atom::i())
            .into(),
    ));
    // println!("permutation{}", gamma_key.rep_permutation);
    weyl.insert_explicit(gamma_conj_key);
    let gamma_adj_key = PermutedStructure::identity(ParamOrConcrete::param(
        gamma_adj_data_weyl(AGS.gamma_adj_strct::<Aind>(4), one.clone(), zero.clone())
            .map_data(|a| a.re + a.im * Atom::i())
            .into(),
    ));
    // println!("permutation{}", gamma_key.rep_permutation);
    weyl.insert_explicit(gamma_adj_key);
    let gamma0_key = PermutedStructure::identity(ParamOrConcrete::param(
        gamma0_weyl(AGS.gamma0_strct::<Aind>(4), one.clone(), zero.clone())
            .map_data(|a| a.re + a.im * Atom::i())
            .into(),
    ));
    // println!("permutation{}", gamma_key.rep_permutation);
    weyl.insert_explicit(gamma0_key);

    let gamma5_key = PermutedStructure::identity(ParamOrConcrete::param(
        gamma5_weyl_data(AGS.gamma5_strct::<Aind>(4), one.clone(), zero.clone())
            .map_data(|a| a.re + a.im * Atom::i())
            .into(),
    ));
    weyl.insert_explicit(gamma5_key);

    let projm_key = PermutedStructure::identity(ParamOrConcrete::param(
        proj_m_data_weyl(AGS.projm_strct::<Aind>(4), one.clone(), zero.clone())
            .map_data(|a| a.re + a.im * Atom::i())
            .into(),
    ));
    weyl.insert_explicit(projm_key);

    let projp_key = PermutedStructure::identity(ParamOrConcrete::param(
        proj_p_data_weyl(AGS.projp_strct::<Aind>(4), one.clone(), zero.clone())
            .map_data(|a| a.re + a.im * Atom::i())
            .into(),
    ));
    weyl.insert_explicit(projp_key);

    weyl
}

pub type HepTensor<Aind> = MixedTensor<f64, ShadowedStructure<Aind>>;

pub type HepNet<Aind> =
    Network<NetworkStore<HepTensor<Aind>, Atom>, ExplicitKey<Aind>, Symbol, Aind>;

pub static HEP_LIB: LazyLock<
    TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>,
> = LazyLock::new(|| hep_lib(1., 0.));

pub static FUN_LIB: LazyLock<
    SymbolLib<RealOrComplexTensor<f64, ShadowedStructure<AbstractIndex>>, PanicMissingConcrete>,
> = LazyLock::new(|| {
    let mut lib = PanicMissingConcrete::new_lib();
    lib.insert(INBUILTS.conj, |a| match a {
        RealOrComplexTensor::Complex(c) => RealOrComplexTensor::Complex(c.map_data(|x| x.conj())),
        RealOrComplexTensor::Real(r) => RealOrComplexTensor::Real(r),
    });
    lib
});

#[cfg(test)]
mod tests {

    use spenso::{
        network::{
            Network, SingleSmallestDegree, SmallestDegreeIter, Steps,
            parsing::{ParseSettings, ShadowedStructure},
            store::NetworkStore,
        },
        structure::{HasStructure, abstract_index::AbstractIndex},
    };
    use symbolica::{
        atom::{Atom, Symbol},
        parse, parse_lit,
    };

    use super::*;

    #[test]
    fn simple_scalar() {
        initialize();
        let _a = HEP_LIB.get(&AGS.gamma_strct(4)).unwrap();

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
        >::try_from_view(expr.as_view(), &*HEP_LIB, &ParseSettings::default())
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
        let _a = HEP_LIB.get(&AGS.gamma_strct(4)).unwrap();

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
        >::try_from_view(expr.as_view(), &*HEP_LIB, &ParseSettings::default())
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
