// use std::{
//     ops::Neg,
//     sync::{LazyLock, RwLock},
// };

// use once_cell::sync::Lazy;
// use spenso::{
//     complex::Complex,
//     data::{SetTensorData, SparseTensor},
//     network::tensor_library::symbolic::{ExplicitKey, TensorLibrary},
//     parametric::MixedTensor,
//     structure::{
//         representation::{LibraryRep, Minkowski, RepName},
//         TensorStructure,
//     },
// };
// use spenso_macros::SimpleRepresentation;
// // use symbolica::{atom::Symbol, symbol};

// #[rustfmt::skip]
// #[derive(SimpleRepresentation)]
// #[derive(
//     Debug,
//     Clone,
//     Copy,
//     PartialEq,
//     Eq,
//     Hash,
//     PartialOrd,
//     Ord,
//     Default,
// )]
// #[representation(name = "bis", self_dual)] // Specify the dual name
// pub struct Bispinor {}

// pub struct GammaLibrary {
//     pub gamma: Symbol,
//     pub projp: Symbol,
//     pub projm: Symbol,
//     pub gamma5: Symbol,
//     pub sigma: Symbol,
// }

// pub static WEYL: LazyLock<GammaLibrary> = LazyLock::new(|| GammaLibrary {
//     gamma: symbol!("weyl::gamma"),
//     projp: symbol!("weyl::projp"),
//     projm: symbol!("weyl::projm"),
//     gamma5: symbol!("weyl::gamma5"),
//     sigma: symbol!("weyl::sigma"),
// });

// pub static WEYLIB: Lazy<RwLock<TensorLibrary<MixedTensor<f64, ExplicitKey>>>> = Lazy::new(|| {
//     let mut library = TensorLibrary::new();

//     let gamma = ExplicitKey::from_iter(
//         [
//             LibraryRep::from(Minkowski {}).new_rep(4),
//             Bispinor {}.new_rep(4).cast(),
//             Bispinor {}.new_rep(4).cast(),
//         ],
//         WEYL.gamma,
//         None,
//     );

//     library.insert_explicit(gamma_data_weyl(gamma, 1., 0.).into());

//     let gamma5 = ExplicitKey::from_iter(
//         [Bispinor {}.new_rep(4), Bispinor {}.new_rep(4)],
//         WEYL.gamma5,
//         None,
//     );
//     library.insert_explicit(gamma5_weyl_data(gamma5, 1., 0.).into());

//     let projm_key = ExplicitKey::from_iter(
//         [Bispinor {}.new_rep(4), Bispinor {}.new_rep(4)],
//         WEYL.projm,
//         None,
//     );
//     library.insert_explicit(proj_m_data_weyl(projm_key, 1., 0.).into());

//     let projp_key = ExplicitKey::from_iter(
//         [Bispinor {}.new_rep(4), Bispinor {}.new_rep(4)],
//         WEYL.projp,
//         None,
//     );
//     library.insert_explicit(proj_p_data_weyl(projp_key, 1., 0.).into());

//     RwLock::new(library)
// });

// #[allow(clippy::similar_names)]
// pub fn gamma_data_dirac<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Clone + Neg<Output = T>,
//     N: TensorStructure,
// {
//     let c1 = Complex::<T>::new(one.clone(), zero.clone());
//     let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
//     let ci = Complex::<T>::new(zero.clone(), one.clone());
//     let cni = Complex::<T>::new(zero.clone(), -one.clone());
//     let mut gamma = SparseTensor::empty(structure);
//     // ! No check on actual structure, should expext mink,bis,bis

//     // dirac gamma matrices

//     gamma.set(&[0, 0, 0], c1.clone()).unwrap();
//     gamma.set(&[0, 1, 1], c1.clone()).unwrap();
//     gamma.set(&[0, 2, 2], cn1.clone()).unwrap();
//     gamma.set(&[0, 3, 3], cn1.clone()).unwrap();

//     gamma.set(&[1, 0, 3], c1.clone()).unwrap();
//     gamma.set(&[1, 1, 2], c1.clone()).unwrap();
//     gamma.set(&[1, 2, 1], cn1.clone()).unwrap();
//     gamma.set(&[1, 3, 0], cn1.clone()).unwrap();

//     gamma.set(&[2, 0, 3], cni.clone()).unwrap();
//     gamma.set(&[2, 1, 2], ci.clone()).unwrap();
//     gamma.set(&[2, 2, 1], ci.clone()).unwrap();
//     gamma.set(&[2, 3, 0], cni.clone()).unwrap();

//     gamma.set(&[3, 0, 2], c1.clone()).unwrap();
//     gamma.set(&[3, 1, 3], cn1.clone()).unwrap();
//     gamma.set(&[3, 2, 0], cn1.clone()).unwrap();
//     gamma.set(&[3, 3, 1], c1.clone()).unwrap();

//     gamma //.to_dense()
// }

// #[allow(clippy::similar_names)]
// pub fn gamma_data_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Neg<Output = T> + Clone,
//     N: TensorStructure,
// {
//     let c1 = Complex::<T>::new(one.clone(), zero.clone());
//     let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
//     let ci = Complex::<T>::new(zero.clone(), one.clone());
//     let cni = Complex::<T>::new(zero.clone(), -one.clone());
//     let mut gamma = SparseTensor::empty(structure);
//     // ! No check on actual structure, should expext mink,bis,bis

//     // dirac gamma matrices

//     gamma.set(&[0, 0, 2], c1.clone()).unwrap();
//     gamma.set(&[0, 1, 3], c1.clone()).unwrap();
//     gamma.set(&[0, 2, 0], c1.clone()).unwrap();
//     gamma.set(&[0, 3, 1], c1.clone()).unwrap();

//     gamma.set(&[1, 0, 3], c1.clone()).unwrap();
//     gamma.set(&[1, 1, 2], c1.clone()).unwrap();
//     gamma.set(&[1, 2, 1], cn1.clone()).unwrap();
//     gamma.set(&[1, 3, 0], cn1.clone()).unwrap();

//     gamma.set(&[2, 0, 3], cni.clone()).unwrap();
//     gamma.set(&[2, 1, 2], ci.clone()).unwrap();
//     gamma.set(&[2, 2, 1], ci.clone()).unwrap();
//     gamma.set(&[2, 3, 0], cni.clone()).unwrap();

//     gamma.set(&[3, 0, 2], c1.clone()).unwrap();
//     gamma.set(&[3, 1, 3], cn1.clone()).unwrap();
//     gamma.set(&[3, 2, 0], cn1.clone()).unwrap();
//     gamma.set(&[3, 3, 1], c1.clone()).unwrap();

//     gamma //.to_dense()
// }

// pub fn gamma0_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Clone,
//     N: TensorStructure,
// {
//     let c1 = Complex::<T>::new(one, zero);
//     let mut gamma0 = SparseTensor::empty(structure);
//     // ! No check on actual structure, should expext bis,bis,lor

//     // dirac gamma0 matrices

//     gamma0.set(&[0, 2], c1.clone()).unwrap();
//     gamma0.set(&[1, 3], c1.clone()).unwrap();
//     gamma0.set(&[2, 0], c1.clone()).unwrap();
//     gamma0.set(&[3, 1], c1.clone()).unwrap();

//     gamma0
// }

// pub fn gamma5_dirac_data<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Clone,
//     N: TensorStructure,
// {
//     let c1 = Complex::<T>::new(one, zero);

//     let mut gamma5 = SparseTensor::empty(structure);

//     gamma5.set(&[0, 2], c1.clone()).unwrap();
//     gamma5.set(&[1, 3], c1.clone()).unwrap();
//     gamma5.set(&[2, 0], c1.clone()).unwrap();
//     gamma5.set(&[3, 1], c1.clone()).unwrap();

//     gamma5
// }

// pub fn gamma5_weyl_data<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Clone + Neg<Output = T>,
//     N: TensorStructure,
// {
//     let c1 = Complex::<T>::new(one, zero);

//     let mut gamma5 = SparseTensor::empty(structure);

//     gamma5.set(&[0, 0], -c1.clone()).unwrap();
//     gamma5.set(&[1, 1], -c1.clone()).unwrap();
//     gamma5.set(&[2, 2], c1.clone()).unwrap();
//     gamma5.set(&[3, 3], c1.clone()).unwrap();

//     gamma5
// }

// #[allow(clippy::similar_names)]
// pub fn proj_m_data_dirac<T, N>(structure: N, half: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Clone + Neg<Output = T>,
//     N: TensorStructure,
// {
//     // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
//     let chalf = Complex::<T>::new(half.clone(), zero.clone());
//     let cnhalf = Complex::<T>::new(-half, zero);

//     let mut proj_m = SparseTensor::empty(structure);

//     proj_m.set(&[0, 0], chalf.clone()).unwrap();
//     proj_m.set(&[1, 1], chalf.clone()).unwrap();
//     proj_m.set(&[2, 2], chalf.clone()).unwrap();
//     proj_m.set(&[3, 3], chalf.clone()).unwrap();

//     proj_m.set(&[0, 2], cnhalf.clone()).unwrap();
//     proj_m.set(&[1, 3], cnhalf.clone()).unwrap();
//     proj_m.set(&[2, 0], cnhalf.clone()).unwrap();
//     proj_m.set(&[3, 1], cnhalf.clone()).unwrap();

//     proj_m
// }

// #[allow(clippy::similar_names)]
// pub fn proj_m_data_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Clone,
//     N: TensorStructure,
// {
//     // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
//     let c1 = Complex::<T>::new(one, zero);
//     let mut proj_m = SparseTensor::empty(structure);

//     proj_m.set(&[0, 0], c1.clone()).unwrap();
//     proj_m.set(&[1, 1], c1.clone()).unwrap();

//     proj_m
// }

// pub fn proj_p_data_dirac<T, N>(structure: N, half: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Clone,
//     N: TensorStructure,
// {
//     // ProjP(1,2) Right chirality projector (( 1+γ5)/ 2 )_s1_s2
//     let chalf = Complex::<T>::new(half, zero);

//     let mut proj_p = SparseTensor::empty(structure);

//     proj_p
//         .set(&[0, 0], chalf.clone())
//         .unwrap_or_else(|_| unreachable!());
//     proj_p
//         .set(&[1, 1], chalf.clone())
//         .unwrap_or_else(|_| unreachable!());
//     proj_p
//         .set(&[2, 2], chalf.clone())
//         .unwrap_or_else(|_| unreachable!());
//     proj_p
//         .set(&[3, 3], chalf.clone())
//         .unwrap_or_else(|_| unreachable!());

//     proj_p
//         .set(&[0, 2], chalf.clone())
//         .unwrap_or_else(|_| unreachable!());
//     proj_p
//         .set(&[1, 3], chalf.clone())
//         .unwrap_or_else(|_| unreachable!());
//     proj_p
//         .set(&[2, 0], chalf.clone())
//         .unwrap_or_else(|_| unreachable!());
//     proj_p
//         .set(&[3, 1], chalf.clone())
//         .unwrap_or_else(|_| unreachable!());

//     proj_p
// }

// #[allow(clippy::similar_names)]
// pub fn proj_p_data_weyl<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Clone,
//     N: TensorStructure,
// {
//     // ProjM(1,2) Left chirality projector (( 1−γ5)/ 2 )_s1_s2
//     let c1 = Complex::<T>::new(one, zero);
//     let mut proj_m = SparseTensor::empty(structure);

//     proj_m.set(&[2, 2], c1.clone()).unwrap();
//     proj_m.set(&[3, 3], c1.clone()).unwrap();

//     proj_m
// }

// #[allow(clippy::similar_names)]
// pub fn sigma_data<T, N>(structure: N, one: T, zero: T) -> SparseTensor<Complex<T>, N>
// where
//     T: Clone + Neg<Output = T>,
//     N: TensorStructure,
// {
//     let c1 = Complex::<T>::new(one.clone(), zero.clone());
//     let cn1 = Complex::<T>::new(-one.clone(), zero.clone());
//     let ci = Complex::<T>::new(zero.clone(), one.clone());
//     let cni = Complex::<T>::new(zero.clone(), -one.clone());

//     let mut sigma = SparseTensor::empty(structure);
//     sigma.set(&[0, 2, 0, 1], c1.clone()).unwrap();
//     sigma.set(&[0, 2, 3, 0], c1.clone()).unwrap();
//     sigma.set(&[0, 3, 1, 2], c1.clone()).unwrap();
//     sigma.set(&[1, 0, 2, 2], c1.clone()).unwrap();
//     sigma.set(&[1, 1, 1, 2], c1.clone()).unwrap();
//     sigma.set(&[1, 3, 0, 2], c1.clone()).unwrap();
//     sigma.set(&[2, 2, 1, 0], c1.clone()).unwrap();
//     sigma.set(&[2, 2, 2, 1], c1.clone()).unwrap();
//     sigma.set(&[2, 3, 3, 2], c1.clone()).unwrap();
//     sigma.set(&[3, 0, 0, 2], c1.clone()).unwrap();
//     sigma.set(&[3, 3, 2, 2], c1.clone()).unwrap();
//     sigma.set(&[3, 1, 3, 2], c1.clone()).unwrap();
//     sigma.set(&[0, 1, 3, 0], ci.clone()).unwrap();
//     sigma.set(&[0, 3, 1, 1], ci.clone()).unwrap();
//     sigma.set(&[0, 3, 2, 0], ci.clone()).unwrap();
//     sigma.set(&[1, 0, 3, 3], ci.clone()).unwrap();
//     sigma.set(&[1, 1, 0, 3], ci.clone()).unwrap();
//     sigma.set(&[1, 1, 2, 0], ci.clone()).unwrap();
//     sigma.set(&[2, 1, 1, 0], ci.clone()).unwrap();
//     sigma.set(&[2, 3, 0, 0], ci.clone()).unwrap();
//     sigma.set(&[2, 3, 3, 1], ci.clone()).unwrap();
//     sigma.set(&[3, 0, 1, 3], ci.clone()).unwrap();
//     sigma.set(&[3, 1, 0, 0], ci.clone()).unwrap();
//     sigma.set(&[3, 1, 2, 3], ci.clone()).unwrap();
//     sigma.set(&[0, 0, 3, 2], cn1.clone()).unwrap();
//     sigma.set(&[0, 1, 0, 2], cn1.clone()).unwrap();
//     sigma.set(&[0, 2, 1, 3], cn1.clone()).unwrap();
//     sigma.set(&[1, 2, 0, 3], cn1.clone()).unwrap();
//     sigma.set(&[1, 2, 1, 1], cn1.clone()).unwrap();
//     sigma.set(&[1, 2, 2, 0], cn1.clone()).unwrap();
//     sigma.set(&[2, 0, 1, 2], cn1.clone()).unwrap();
//     sigma.set(&[2, 1, 2, 2], cn1.clone()).unwrap();
//     sigma.set(&[2, 2, 3, 3], cn1.clone()).unwrap();
//     sigma.set(&[3, 2, 0, 0], cn1.clone()).unwrap();
//     sigma.set(&[3, 2, 2, 3], cn1.clone()).unwrap();
//     sigma.set(&[3, 2, 3, 1], cn1.clone()).unwrap();
//     sigma.set(&[0, 0, 2, 3], cni.clone()).unwrap();
//     sigma.set(&[0, 0, 3, 1], cni.clone()).unwrap();
//     sigma.set(&[0, 1, 1, 3], cni.clone()).unwrap();
//     sigma.set(&[1, 0, 2, 1], cni.clone()).unwrap();
//     sigma.set(&[1, 3, 0, 1], cni.clone()).unwrap();
//     sigma.set(&[1, 3, 3, 0], cni.clone()).unwrap();
//     sigma.set(&[2, 0, 0, 3], cni.clone()).unwrap();
//     sigma.set(&[2, 0, 1, 1], cni.clone()).unwrap();
//     sigma.set(&[2, 1, 3, 3], cni.clone()).unwrap();
//     sigma.set(&[3, 0, 0, 1], cni.clone()).unwrap();
//     sigma.set(&[3, 3, 1, 0], cni.clone()).unwrap();
//     sigma.set(&[3, 3, 2, 1], cni.clone()).unwrap();

//     sigma
// }
