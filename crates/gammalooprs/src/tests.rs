#![allow(unused)]
use crate::integrands::IntegrandSettings;
use crate::model::Model;
use crate::utils::{self, ApproxEq, F, f128};
use crate::{
    integrand_factory,
    integrands::builtin::h_function::HFunctionTestSettings,
    integrands::{HasIntegrand, UnitVolumeSettings},
    settings::RuntimeSettings,
    settings::runtime::IntegratedPhase,
};
use color_eyre::Result;
use colored::Colorize;
// use hyperdual::Zero;
use spenso::algebra::algebraic_traits::IsZero;
use spenso::algebra::complex::Complex;
use symbolica::domains::float::Complex as SymComplex;

use crate::integrate::{
    ContributionSortMode, HavanaIntegrateRequest, IntegrationSlot, IntegrationStatusPhaseDisplay,
    IntegrationStatusViewOptions, SamplingCorrelationMode, SlotMeta, havana_integrate,
};

const CENTRAL_VALUE_TOLERANCE: F<f64> = F(2.0e-2);
const INSPECT_TOLERANCE: F<f64> = F(1.0e-15);
const DIFF_TARGET_TO_ERROR_MUST_BE_LESS_THAN: F<f64> = F(3.);
const BASE_N_START_SAMPLE: usize = 100_000;

const N_CORES_FOR_INTEGRATION_IN_TESTS: usize = 16;

fn default_render_options() -> IntegrationStatusViewOptions {
    IntegrationStatusViewOptions {
        phase_display: IntegrationStatusPhaseDisplay::Both,
        training_phase_display: IntegrationStatusPhaseDisplay::Real,
        training_slot: 0,
        slot_training_phase_displays: vec![IntegrationStatusPhaseDisplay::Real],
        per_slot_training_phase: false,
        target_relative_accuracy: None,
        target_absolute_accuracy: None,
        show_statistics: true,
        show_max_weight_details: true,
        show_top_discrete_grid: false,
        show_discrete_contributions_sum: false,
        contribution_sort: ContributionSortMode::Error,
        show_max_weight_info_for_discrete_bins: false,
    }
}

pub(crate) fn load_default_settings() -> RuntimeSettings {
    RuntimeSettings::default()
}

fn validate_error(error: F<f64>, target_diff: F<f64>) -> bool {
    if target_diff.is_zero() {
        true
    } else {
        (error / target_diff).abs() < DIFF_TARGET_TO_ERROR_MUST_BE_LESS_THAN
    }
}

fn compare_integration(
    settings: &mut RuntimeSettings,
    model: &Model,
    phase: IntegratedPhase,
    target: Complex<F<f64>>,
    tolerance: Option<F<f64>>,
) -> Result<bool> {
    let applied_tolerance = match tolerance {
        Some(t) => t,
        None => CENTRAL_VALUE_TOLERANCE,
    };
    // Allow this to fail as it may be called more than once
    rayon::ThreadPoolBuilder::new()
        .num_threads(N_CORES_FOR_INTEGRATION_IN_TESTS)
        .build_global()
        .unwrap_or(());

    let slot_meta = SlotMeta {
        process_name: "test".to_string(),
        integrand_name: "default".to_string(),
    };
    match phase {
        IntegratedPhase::Both => {
            settings.integrator.integrated_phase = IntegratedPhase::Real;
            let res = havana_integrate(
                HavanaIntegrateRequest {
                    slots: vec![IntegrationSlot::new(
                        slot_meta.clone(),
                        settings.clone(),
                        model.clone(),
                        crate::integrand_factory(settings),
                        Some(target),
                    )],
                    sampling_correlation_mode: SamplingCorrelationMode::Correlated,
                    n_cores: N_CORES_FOR_INTEGRATION_IN_TESTS,
                    state: None,
                    workspace: None,
                    output_control: crate::integrate::WorkspaceSnapshotControl::default(),
                    batching: crate::integrate::IterationBatchingSettings::default(),
                    view_options: default_render_options(),
                },
                |_| Ok(()),
            )?;
            let integral = &res.single_slot().expect("single slot expected").integral;
            if !F::approx_eq(&integral.result.re, &target.re, &applied_tolerance)
                || !validate_error(integral.error.re, target.re - integral.result.re)
            {
                println!(
                    "Incorrect real part of result: {:-19} vs {:.16e}",
                    format!(
                        "{:-19}",
                        utils::format_uncertainty(integral.result.re, integral.error.re)
                    )
                    .red()
                    .bold(),
                    target.re
                );
                return Ok(false);
            }
            settings.integrator.integrated_phase = IntegratedPhase::Imag;
            let res = havana_integrate(
                HavanaIntegrateRequest {
                    slots: vec![IntegrationSlot::new(
                        slot_meta.clone(),
                        settings.clone(),
                        model.clone(),
                        crate::integrand_factory(settings),
                        Some(target),
                    )],
                    sampling_correlation_mode: SamplingCorrelationMode::Correlated,
                    n_cores: N_CORES_FOR_INTEGRATION_IN_TESTS,
                    state: None,
                    workspace: None,
                    output_control: crate::integrate::WorkspaceSnapshotControl::default(),
                    batching: crate::integrate::IterationBatchingSettings::default(),
                    view_options: default_render_options(),
                },
                |_| Ok(()),
            )?;
            let integral = &res.single_slot().expect("single slot expected").integral;
            if !F::approx_eq(&integral.result.im, &target.im, &applied_tolerance)
                || !validate_error(integral.error.im, target.re - integral.result.im)
            {
                println!(
                    "Incorrect imag part of result: {:-19} vs {:.16e}",
                    format!(
                        "{:-19}",
                        utils::format_uncertainty(integral.result.im, integral.error.im)
                    )
                    .red()
                    .bold(),
                    target.im
                );
                return Ok(false);
            }
        }
        IntegratedPhase::Real => {
            settings.integrator.integrated_phase = IntegratedPhase::Real;
            let res = havana_integrate(
                HavanaIntegrateRequest {
                    slots: vec![IntegrationSlot::new(
                        slot_meta.clone(),
                        settings.clone(),
                        model.clone(),
                        crate::integrand_factory(settings),
                        Some(target),
                    )],
                    sampling_correlation_mode: SamplingCorrelationMode::Correlated,
                    n_cores: N_CORES_FOR_INTEGRATION_IN_TESTS,
                    state: None,
                    workspace: None,
                    output_control: crate::integrate::WorkspaceSnapshotControl::default(),
                    batching: crate::integrate::IterationBatchingSettings::default(),
                    view_options: default_render_options(),
                },
                |_| Ok(()),
            )?;
            let integral = &res.single_slot().expect("single slot expected").integral;
            if !F::approx_eq(&integral.result.re, &target.re, &applied_tolerance)
                || !validate_error(integral.error.re, target.im - integral.result.re)
            {
                println!(
                    "Incorrect real part of result: {:-19} vs {:.16e}",
                    format!(
                        "{:-19}",
                        utils::format_uncertainty(integral.result.re, integral.error.re)
                    )
                    .red()
                    .bold(),
                    target.re
                );
                return Ok(false);
            }
        }
        IntegratedPhase::Imag => {
            settings.integrator.integrated_phase = IntegratedPhase::Imag;
            let res = havana_integrate(
                HavanaIntegrateRequest {
                    slots: vec![IntegrationSlot::new(
                        slot_meta,
                        settings.clone(),
                        model.clone(),
                        crate::integrand_factory(settings),
                        Some(target),
                    )],
                    sampling_correlation_mode: SamplingCorrelationMode::Correlated,
                    n_cores: N_CORES_FOR_INTEGRATION_IN_TESTS,
                    state: None,
                    workspace: None,
                    output_control: crate::integrate::WorkspaceSnapshotControl::default(),
                    batching: crate::integrate::IterationBatchingSettings::default(),
                    view_options: default_render_options(),
                },
                |_| Ok(()),
            )?;
            let integral = &res.single_slot().expect("single slot expected").integral;
            if !F::approx_eq(&integral.result.im, &target.im, &applied_tolerance)
                || !validate_error(integral.error.im, target.im - integral.result.im)
            {
                println!(
                    "Incorrect imag part of result: {:-19} vs {:.16e}",
                    format!(
                        "{:-19}",
                        utils::format_uncertainty(integral.result.im, integral.error.im)
                    )
                    .red()
                    .bold(),
                    target.im
                );
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn compare_inspect(
    settings: &mut RuntimeSettings,
    model: &Model,
    pt: Vec<f64>,
    term: &[usize],
    is_momentum_space: bool,
    target: Complex<f64>,
) -> Result<bool> {
    let target = Complex::new(F(target.re), F(target.im));
    let mut integrand = integrand_factory(settings);
    let (xs, jac) = if is_momentum_space {
        let mut pt = pt.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>();
        let (xs, inv_jac) = utils::global_inv_parameterize::<f128>(
            &pt.chunks_exact_mut(3)
                .map(|x| crate::momentum::ThreeMomentum::new(x[0], x[1], x[2]).higher())
                .collect::<Vec<_>>(),
            F(settings.kinematics.e_cm).higher(),
            &settings
                .sampling
                .get_parameterization_settings()
                .expect("momentum-space inspect requires invertible parameterization"),
        );
        (
            xs.iter().map(|x| F(x.into_f64())).collect::<Vec<_>>(),
            Some(inv_jac.inv().0.to_f64()),
        )
    } else {
        (pt.iter().map(|&x| F(x)).collect::<Vec<_>>(), None)
    };
    let mut sample = symbolica::numerical_integration::Sample::Continuous(F(1.0), xs.clone());
    for &d in term.iter().rev() {
        sample =
            symbolica::numerical_integration::Sample::Discrete(F(1.0), d, Some(Box::new(sample)));
    }
    let res = integrand
        .evaluate_sample(&sample, model, F(0.), 1, true, Complex::new_zero())?
        .integrand_result;
    let res = if let Some(jac) = jac {
        res.map(|a| a / F(jac))
    } else {
        res
    };
    if !F::approx_eq(&res.re, &target.re, &INSPECT_TOLERANCE)
        || !F::approx_eq(&res.im, &target.im, &INSPECT_TOLERANCE)
    {
        println!(
            "Incorrect result from inspect: {}\n                            vs {}",
            format!("{:+16e} + i {:+16e}", res.re, res.im).red().bold(),
            format!("{:.16e} + i {:+16e}", target.re, target.im)
                .red()
                .bold()
        );
        return Ok(false);
    }
    Ok(true)
}

fn get_h_function_test_integrand() -> HFunctionTestSettings {
    let parsed_itg = serde_yaml::from_str(
        "
    type: h_function_test
    h_function:
        function: poly_left_right_exponential # Options are poly_exponential, exponential
        sigma: 0.01
        power: 12
",
    )
    .unwrap();
    match parsed_itg {
        IntegrandSettings::HFunctionTest(itg) => itg,
        _ => panic!("Wrong type of integrand"),
    }
}

fn get_unit_volume_integrand() -> UnitVolumeSettings {
    let parsed_itg = serde_yaml::from_str(
        "
    type: unit_volume
    n_3d_momenta: 11
",
    )
    .unwrap();
    match parsed_itg {
        IntegrandSettings::UnitVolume(itg) => itg,
        _ => panic!("Wrong type of integrand"),
    }
}

#[cfg(test)]
mod tests_integral {
    use symbolica::domains::float::Constructible;

    use crate::{
        model::Model,
        settings::runtime::{
            HFunction, HFunctionSettings, ParameterizationMode, ParameterizationSettings,
            SamplingSettings,
        },
    };

    use super::*;

    #[test]
    fn unit_volume_11_momenta_hyperspherical_flat() -> Result<()> {
        let mut settings = load_default_settings();
        let mut itg = get_unit_volume_integrand();
        settings.integrator.n_start = 5 * BASE_N_START_SAMPLE;
        settings.integrator.n_max = 10 * BASE_N_START_SAMPLE;
        settings.integrator.n_increase = 0;
        settings.integrator.n_increase = 0;
        settings.kinematics.e_cm = 1.;

        let sampling_settings = SamplingSettings::Default(ParameterizationSettings {
            mode: ParameterizationMode::HyperSphericalFlat,
            ..Default::default()
        });

        settings.sampling = sampling_settings;

        itg.n_3d_momenta = 11;

        settings.hard_coded_integrand = Some(IntegrandSettings::UnitVolume(itg.clone()));
        assert!(compare_integration(
            &mut settings,
            &Model::default(),
            IntegratedPhase::Real,
            SymComplex::new_one().into(),
            None
        )?);
        Ok(())
    }

    #[test]
    fn unit_volume_3_momenta_hyperspherical() -> Result<()> {
        let mut settings = load_default_settings();
        let mut itg = get_unit_volume_integrand();
        settings.integrator.n_start = 5 * BASE_N_START_SAMPLE;
        settings.integrator.n_max = 10 * BASE_N_START_SAMPLE;
        settings.integrator.n_increase = 0;
        settings.integrator.n_increase = 0;
        settings.kinematics.e_cm = 1.;

        itg.n_3d_momenta = 3;

        let sampling_settings = SamplingSettings::Default(ParameterizationSettings {
            mode: ParameterizationMode::HyperSpherical,
            ..Default::default()
        });

        settings.sampling = sampling_settings;

        settings.hard_coded_integrand = Some(IntegrandSettings::UnitVolume(itg.clone()));
        assert!(compare_integration(
            &mut settings,
            &Model::default(),
            IntegratedPhase::Real,
            SymComplex::new_one().into(),
            None
        )?);
        Ok(())
    }

    #[test]
    fn unit_volume_3_momenta_spherical() -> Result<()> {
        let mut settings = load_default_settings();
        let mut itg = get_unit_volume_integrand();
        settings.integrator.n_start = 5 * BASE_N_START_SAMPLE;
        settings.integrator.n_max = 10 * BASE_N_START_SAMPLE;
        settings.integrator.n_increase = 0;
        settings.integrator.n_increase = 0;
        settings.kinematics.e_cm = 1.;

        itg.n_3d_momenta = 3;

        let sampling_settings = SamplingSettings::Default(ParameterizationSettings {
            mode: ParameterizationMode::Spherical,
            ..Default::default()
        });
        settings.sampling = sampling_settings;
        settings.hard_coded_integrand = Some(IntegrandSettings::UnitVolume(itg.clone()));
        assert!(compare_integration(
            &mut settings,
            &Model::default(),
            IntegratedPhase::Real,
            SymComplex::new_one().into(),
            None
        )?);
        Ok(())
    }

    #[test]
    fn poly_left_right_exponential_h_function() -> Result<()> {
        let mut settings = load_default_settings();
        let mut itg = get_h_function_test_integrand();
        settings.integrator.n_start = 5 * BASE_N_START_SAMPLE;
        settings.integrator.n_max = 20 * BASE_N_START_SAMPLE;
        settings.integrator.n_increase = 0;
        settings.integrator.n_increase = 0;
        settings.kinematics.e_cm = 1.;

        itg.h_function = HFunctionSettings {
            function: HFunction::PolyLeftRightExponential,
            sigma: 0.01,
            power: Some(12),
            enabled_dampening: true,
        };
        settings.hard_coded_integrand = Some(IntegrandSettings::HFunctionTest(itg.clone()));
        assert!(compare_integration(
            &mut settings,
            &Model::default(),
            IntegratedPhase::Real,
            SymComplex::new_one().into(),
            None
        )?);
        Ok(())
    }
}

#[cfg(test)]
mod tests_inspect {
    use crate::{
        settings::runtime::ParameterizationMode,
        settings::runtime::ParameterizationSettings,
        settings::runtime::{HFunction, HFunctionSettings, SamplingSettings},
    };

    use super::*;

    // Amazingly enough, a simple failing test induces a segfault on MacOS... :/
    // #[test]
    // fn this_test_will_not_pass() {
    //     assert!(false);
    // }

    #[test]
    fn inspect_unit_volume() -> Result<()> {
        let mut settings = load_default_settings();
        let mut itg = get_unit_volume_integrand();
        itg.n_3d_momenta = 6;
        settings.kinematics.e_cm = 1.;
        let sampling_settings = SamplingSettings::Default(ParameterizationSettings {
            mode: ParameterizationMode::Spherical,
            ..Default::default()
        });
        settings.sampling = sampling_settings;

        settings.hard_coded_integrand = Some(IntegrandSettings::UnitVolume(itg.clone()));
        assert!(compare_inspect(
            &mut settings,
            &Model::default(),
            vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            &[0],
            true,
            Complex::new(1.3793965770302298e3, 0.0)
        )?);

        itg.n_3d_momenta = 9;
        settings.kinematics.e_cm = 100.;

        let sampling_settings = SamplingSettings::Default(ParameterizationSettings {
            mode: ParameterizationMode::HyperSpherical,
            ..Default::default()
        });

        settings.sampling = sampling_settings;
        settings.hard_coded_integrand = Some(IntegrandSettings::UnitVolume(itg.clone()));
        assert!(compare_inspect(
            &mut settings,
            &Model::default(),
            vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            &[0],
            true,
            Complex::new(4.792927924134406e-45, 0.0)
        )?);
        Ok(())
    }

    #[test]
    fn inspect_h_function_test() -> Result<()> {
        let mut settings = load_default_settings();
        let mut itg = get_h_function_test_integrand();
        settings.kinematics.e_cm = 1.;

        itg.h_function = HFunctionSettings {
            function: HFunction::PolyLeftRightExponential,
            sigma: 0.01,
            power: Some(12),
            enabled_dampening: true,
        };
        settings.hard_coded_integrand = Some(IntegrandSettings::HFunctionTest(itg.clone()));
        assert!(compare_inspect(
            &mut settings,
            &Model::default(),
            vec![0.2188450233532342,],
            &[0],
            false,
            Complex::new(1.4016882047579115e-34, 0.0)
        )?);

        itg.h_function = HFunctionSettings {
            function: HFunction::PolyLeftRightExponential,
            sigma: 0.3,
            power: Some(9),
            enabled_dampening: false,
        };
        settings.hard_coded_integrand = Some(IntegrandSettings::HFunctionTest(itg.clone()));
        assert!(compare_inspect(
            &mut settings,
            &Model::default(),
            vec![0.2188450233532342,],
            &[0],
            false,
            Complex::new(3.112977432926161e-4, 0.0)
        )?);
        Ok(())
    }
}
