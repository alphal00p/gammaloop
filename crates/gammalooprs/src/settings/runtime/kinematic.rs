pub mod improvement;

use bincode_trait_derive::{Decode, Encode};
use eyre::eyre;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use tabled::{builder::Builder, settings::Style};
use tracing::debug;
use typed_index_collections::TiVec;

use crate::{
    DependentMomentaConstructor, GammaLoopContext,
    graph::Graph,
    momentum::sample::{ExternalFourMomenta, ExternalIndex},
    momentum::signature::ExternalSignature,
    momentum::{
        self, Dep, ExternalMomenta, FourMomentum, Helicity, Polarization, Rotatable, SignOrZero,
    },
    settings::runtime::kinematic::improvement::{PhaseSpaceImprovementSettings, improve_ps},
    utils::{
        F, FloatLike, f128,
        serde_utils::{IsDefault, is_float},
    },
};
use color_eyre::{Result, Section};

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, PartialEq, JsonSchema)]
#[trait_decode(trait= GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
pub struct KinematicsSettings {
    #[serde(skip_serializing_if = "is_float::<64>")]
    pub e_cm: f64,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub externals: Externals,
}

impl KinematicsSettings {
    pub fn random(graph: &Graph, seed: u64) -> Self {
        Self {
            e_cm: 64.,
            externals: graph.random_externals(seed),
        }
    }
}

impl Default for KinematicsSettings {
    fn default() -> Self {
        Self {
            e_cm: 64.,
            externals: Externals::default(),
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Encode, Decode, JsonSchema)]
// #[trait_decode(trait= GammaLoopContext)]
#[serde(tag = "type", content = "data")]
#[serde(deny_unknown_fields)]
pub enum Externals {
    #[serde(rename = "constant")]
    Constant {
        momenta: Vec<ExternalMomenta<F<f64>>>,
        helicities: Vec<Helicity>,
        #[serde(default, skip_serializing_if = "IsDefault::is_default")]
        improvement_settings: PhaseSpaceImprovementSettings,
        // must be set in warmup, but is tied to the constant
        #[serde(skip)]
        f_64_cache: Option<TiVec<ExternalIndex, FourMomentum<F<f64>>>>,
        #[serde(skip)]
        f_128_cache: Option<TiVec<ExternalIndex, FourMomentum<F<f128>>>>,
    },
    // add different type of pdfs here when needed
}

impl Rotatable for Externals {
    fn rotate(&self, rotation: &momentum::Rotation) -> Self {
        match self {
            Externals::Constant {
                momenta,
                helicities,
                f_64_cache,
                f_128_cache,
                improvement_settings,
            } => {
                let momenta = momenta.iter().map(|m| m.rotate(rotation)).collect();
                Externals::Constant {
                    momenta,
                    helicities: helicities.clone(),
                    improvement_settings: improvement_settings.clone(),
                    f_64_cache: f_64_cache
                        .as_ref()
                        .map(|cache| cache.iter().map(|m| m.rotate(rotation)).collect()),
                    f_128_cache: f_128_cache
                        .as_ref()
                        .map(|cache| cache.iter().map(|m| m.rotate(rotation)).collect()),
                }
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub enum Polarizations {
    Constant {
        polarizations: Vec<Polarization<Complex<F<f64>>>>,
    },
    None,
}

impl Rotatable for Polarizations {
    fn rotate(&self, rotation: &momentum::Rotation) -> Self {
        match self {
            Polarizations::Constant { polarizations } => {
                let polarizations = polarizations.iter().map(|p| p.rotate(rotation)).collect();
                Polarizations::Constant { polarizations }
            }
            Polarizations::None => Polarizations::None,
        }
    }
}

use thiserror::Error;

#[derive(Error, Debug)]
pub enum ExternalsValidationError {
    #[error("There should be exactly one dependent external momentum")]
    WrongNumberOfDependentMomenta,
    #[error("Found {0} momenta, expected {1}")]
    WrongNumberOfMomentaExpected(usize, usize),
    #[error("Found {0} helicities, expected {1}")]
    WrongNumberOfHelicities(usize, usize),
    #[error("Massless vector cannot have zero helicity: pos {0}")]
    MasslessVectorZeroHelicity(usize),
    #[error("Spinors cannot have zero helicity at pos {0}")]
    SpinorZeroHelicity(usize),
    #[error("Scalars cannot have non-zero helicity at pos {0}")]
    ScalarNonZeroHelicity(usize),
    #[error("{0} is an Unsuported external spin for pos {0}")]
    UnsupportedSpin(isize, usize),
}

impl Externals {
    pub fn validate_helicities(
        &self,
        spins: &[(isize, bool)],
    ) -> Result<(), ExternalsValidationError> {
        match self {
            Externals::Constant { helicities, .. } => {
                if helicities.len() == spins.len() {
                    for (i, (h, (s, is_massless))) in
                        helicities.iter().zip(spins.iter()).enumerate()
                    {
                        match *s {
                            1 => {
                                if !h.is_zero() {
                                    return Err(ExternalsValidationError::ScalarNonZeroHelicity(i));
                                }
                            }
                            2 => {
                                if h.is_zero() {
                                    return Err(ExternalsValidationError::SpinorZeroHelicity(i));
                                }
                            }
                            3 => {
                                if h.is_zero() && *is_massless {
                                    return Err(
                                        ExternalsValidationError::MasslessVectorZeroHelicity(i),
                                    );
                                }
                            }
                            s => return Err(ExternalsValidationError::UnsupportedSpin(s, i)),
                        }
                    }
                    Ok(())
                } else {
                    Err(ExternalsValidationError::WrongNumberOfHelicities(
                        helicities.len(),
                        spins.len(),
                    ))
                }
            }
        }
    }

    pub fn set_dependent_at_end(
        &mut self,
        signature: &ExternalSignature,
    ) -> Result<(), ExternalsValidationError> {
        match self {
            Externals::Constant { momenta, .. } => {
                let mut sum: FourMomentum<F<f128>> = FourMomentum::from([F(0.0); 4]).higher();
                let mut pos_dep = 0;
                let mut n_dep = 0;

                let mut dependent_sign = SignOrZero::Plus;

                for ((i, m), s) in momenta.iter().enumerate().zip(signature.iter()) {
                    if let Ok(a) = FourMomentum::try_from(*m) {
                        sum -= *s * a.higher();
                    } else {
                        pos_dep = i;
                        n_dep += 1;
                        dependent_sign = *s;
                    }
                }
                if n_dep == 1 {
                    momenta[pos_dep] = (dependent_sign * sum.lower()).into();
                    let len = momenta.len();
                    momenta[len - 1] = ExternalMomenta::Dependent(Dep::Dep);
                } else if n_dep == 0 {
                    debug!("No dependent momentum found, adding the sum at the end");
                    momenta.push(ExternalMomenta::Dependent(Dep::Dep));
                } else {
                    return Err(ExternalsValidationError::WrongNumberOfDependentMomenta);
                }

                let len = momenta.len();
                if len == signature.len() {
                    Ok(())
                } else {
                    Err(ExternalsValidationError::WrongNumberOfMomentaExpected(
                        len - 1,
                        signature.len() - 1,
                    ))
                }
            }
        }
    }

    #[inline(never)]
    // #[comemo::memoize]
    pub fn get_dependent_externals<T: FloatLike>(
        &self,
        dependent_momenta_constructor: DependentMomentaConstructor,
    ) -> Result<ExternalFourMomenta<F<T>>>
// where
    //     T::Higher: PrecisionUpgradable<Lower = T> + FloatLike,
    {
        if let Some(cached) = T::try_extract_externals_from_cache(self) {
            return Ok(cached.clone());
        }

        match self {
            Externals::Constant { momenta, .. } => {
                match dependent_momenta_constructor {
                    DependentMomentaConstructor::Amplitude(external_signature) => {
                        if external_signature.is_empty() {
                            return Ok(vec![].into());
                        }

                        let mut sum: FourMomentum<F<T>> = FourMomentum::from([
                            F::<T>::from_f64(0.0),
                            F::from_f64(0.0),
                            F::from_f64(0.0),
                            F::from_f64(0.0),
                        ]);
                        // .higher();
                        let mut pos_dep = external_signature.len() - 1;

                        let mut dependent_sign = SignOrZero::Plus;

                        let mut dependent_momenta = vec![];

                        if momenta.len() != external_signature.len() {
                            return Err(eyre!(
                                "External Momentum in input do not match the number of externals"
                            ))
                            .with_note(|| {
                                let mut table = Builder::new();
                                for m in momenta {
                                    match m {
                                        ExternalMomenta::Dependent(_) => {
                                            table.push_record(["Dependent"])
                                        }
                                        ExternalMomenta::Independent([e, x, y, z]) => table
                                            .push_record([
                                                e.to_string(),
                                                x.to_string(),
                                                y.to_string(),
                                                z.to_string(),
                                            ]),
                                    }
                                }
                                format!(
                                    "External momenta: \n{}\n{}",
                                    table.build().with(Style::rounded()),
                                    external_signature
                                )
                            });
                        }

                        for ((i, m), s) in momenta.iter().enumerate().zip(external_signature.iter())
                        {
                            if let Ok(a) = FourMomentum::try_from(*m) {
                                // println!("external{i}: {}", a);
                                let a = FourMomentum::<F<T>>::from_ff64(&a);
                                sum -= *s * a.clone(); //.higher();
                                dependent_momenta.push(a);
                            } else {
                                pos_dep = i;
                                dependent_sign = *s;
                                dependent_momenta.push(sum.clone()); //.lower());
                            }
                        }

                        dependent_momenta[pos_dep] = dependent_sign * sum; //.lower();
                        let res = dependent_momenta.into();

                        Ok(res)
                    }
                    //DependentMomentaConstructor::CrossSection {
                    //    external_connections,
                    //} => {
                    //    assert_eq!(external_connections.len(), momenta.len());

                    //    // set all_to_zero
                    //    let mut dependent_momenta = (0..momenta.len() * 2)
                    //        .map(|_| {
                    //            FourMomentum::<F<T>>::from_ff64(&FourMomentum::from([F(0.0); 4]))
                    //        })
                    //        .collect::<ExternalFourMomenta<F<T>>>();

                    //    let incoming_momenta = momenta
                    //        .iter()
                    //        .map(|m| {
                    //            FourMomentum::<F<T>>::from_ff64(
                    //                &FourMomentum::try_from(*m)
                    //                    .expect("dependent momenta in cross section not allowed"),
                    //            )
                    //        })
                    //        .collect_vec();

                    //    for (external_connection, incoming_momentum) in
                    //        external_connections.iter().zip(incoming_momenta)
                    //    {
                    //        dependent_momenta[external_connection.incoming_index] =
                    //            incoming_momentum.clone();
                    //        dependent_momenta[external_connection.outgoing_index] =
                    //            incoming_momentum;
                    //    }

                    //    Ok(dependent_momenta)
                    //}
                    DependentMomentaConstructor::CrossSection => {
                        let dependent_momenta = momenta
                            .iter()
                            .map(|m| {
                                FourMomentum::<F<T>>::from_ff64(
                                    &FourMomentum::try_from(*m)
                                        .expect("dependent momenta in None not allowed"),
                                )
                            })
                            .collect::<ExternalFourMomenta<F<T>>>();
                        Ok(dependent_momenta)
                    }
                }
            }
        }
    }

    #[allow(unused_variables)]
    #[inline]
    pub fn get_indep_externals(&self) -> Vec<FourMomentum<F<f64>>> {
        match self {
            Externals::Constant {
                momenta,
                helicities,
                ..
            } => {
                let momenta: Vec<FourMomentum<_>> = momenta
                    .iter()
                    .flat_map(|e| FourMomentum::try_from(*e))
                    .collect();
                momenta
            }
        }
    }

    pub(crate) fn get_helicities(&self) -> &[Helicity] {
        match self {
            Externals::Constant { helicities, .. } => helicities,
        }
    }

    pub(crate) fn _pdf(&self, _x_space_point: &[F<f64>]) -> F<f64> {
        match self {
            Externals::Constant { .. } => F(1.0),
        }
    }

    pub fn improve_and_cache(
        &mut self,
        constructor: DependentMomentaConstructor,
        masses: &TiVec<ExternalIndex, F<f64>>,
        e_cm: &F<f64>,
    ) -> Result<()> {
        let dep_momenta_f64 = self.get_dependent_externals::<f64>(constructor)?;
        let dep_momenta_f128 = self.get_dependent_externals::<f128>(constructor)?;

        match constructor {
            DependentMomentaConstructor::Amplitude(signature) => match self {
                Externals::Constant {
                    improvement_settings,
                    f_64_cache,
                    f_128_cache,
                    ..
                } => {
                    let improved_f64 = improve_ps(
                        &dep_momenta_f64,
                        masses,
                        signature,
                        e_cm,
                        improvement_settings,
                    )?;

                    let upcasted_masses = masses.iter().map(|m| F::<f128>::from_ff64(*m)).collect();

                    let improved_f128 = improve_ps(
                        &dep_momenta_f128,
                        &upcasted_masses,
                        signature,
                        &F::<f128>::from_ff64(*e_cm),
                        improvement_settings,
                    )?;
                    *f_64_cache = Some(improved_f64);
                    *f_128_cache = Some(improved_f128);
                    Ok(())
                }
            },
            // todo, improve the on-shellness
            DependentMomentaConstructor::CrossSection { .. } => Ok(()),
        }
    }
}

#[test]
fn external_inv() {
    let mut ext = Externals::Constant {
        momenta: vec![[F(1.), F(2.), F(3.), F(4.)].into(); 3],
        helicities: vec![Helicity::Plus; 4],
        f_64_cache: None,
        f_128_cache: None,
        improvement_settings: PhaseSpaceImprovementSettings::default(),
    };

    let signs: ExternalSignature = [1i8, 1, 1, 1].into_iter().collect();
    ext.set_dependent_at_end(&signs).unwrap();

    let momenta = vec![
        ExternalMomenta::Dependent(Dep::Dep),
        [F(1.), F(2.), F(3.), F(4.)].into(),
        [F(1.), F(2.), F(3.), F(4.)].into(),
        [F(-3.), F(-6.), F(-9.), F(-12.)].into(),
    ];
    let mut ext2 = Externals::Constant {
        momenta,
        helicities: vec![Helicity::Plus; 4],
        f_64_cache: None,
        f_128_cache: None,
        improvement_settings: PhaseSpaceImprovementSettings::default(),
    };

    ext2.set_dependent_at_end(&signs).unwrap();

    assert_eq!(ext, ext2);
}

impl Default for Externals {
    fn default() -> Self {
        Externals::Constant {
            momenta: vec![],
            helicities: vec![],
            f_64_cache: None,
            f_128_cache: None,
            improvement_settings: PhaseSpaceImprovementSettings::default(),
        }
    }
}
