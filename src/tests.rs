#![allow(dead_code)]
use crate::*;
#[allow(unused)]
use crate::{
    h_function_test::HFunctionTestSettings, integrands::UnitVolumeSettings,
    observables::JetSliceSettings, observables::PhaseSpaceSelectorSettings,
};

use crate::inspect::inspect;
use crate::integrate::{havana_integrate, UserData};

const CENTRAL_VALUE_TOLERANCE: f64 = 2.0e-2;
const INSPECT_TOLERANCE: f64 = 1.0e-15;
const DIFF_TARGET_TO_ERROR_MUST_BE_LESS_THAN: f64 = 3.;
const BASE_N_START_SAMPLE: usize = 100_000;

const N_CORES_FOR_INTEGRATION_IN_TESTS: usize = 16;

fn load_default_settings() -> Settings {
    Settings::from_file("./src/test_resources/default_tests_config.yaml").unwrap()
}

fn approx_eq(res: f64, target: f64, tolerance: f64) -> bool {
    if target == 0.0 {
        return res.abs() < tolerance;
    } else {
        ((res - target) / target).abs() < tolerance
    }
}

fn validate_error(error: f64, target_diff: f64) -> bool {
    if target_diff == 0.0 {
        return true;
    } else {
        (error / target_diff).abs() < DIFF_TARGET_TO_ERROR_MUST_BE_LESS_THAN
    }
}

fn compare_integration(
    settings: &mut Settings,
    phase: IntegratedPhase,
    target: Complex<f64>,
    tolerance: Option<f64>,
) -> bool {
    let applied_tolerance = match tolerance {
        Some(t) => t,
        None => CENTRAL_VALUE_TOLERANCE,
    };
    // Allow this to fail as it may be called more than once
    rayon::ThreadPoolBuilder::new()
        .num_threads(N_CORES_FOR_INTEGRATION_IN_TESTS)
        .build_global()
        .unwrap_or_else(|_| {});

    let user_data_generator = |settings: &Settings| UserData {
        integrand: (0..N_CORES_FOR_INTEGRATION_IN_TESTS)
            .map(|_i| crate::integrand_factory(settings))
            .collect(),
    };
    match phase {
        IntegratedPhase::Both => {
            (*settings).integrator.integrated_phase = IntegratedPhase::Real;
            let res = havana_integrate(&settings, user_data_generator, Some(target));
            if !approx_eq(res.result[0], target.re, applied_tolerance)
                || !validate_error(res.error[0], target.re - res.result[0])
            {
                println!(
                    "Incorrect real part of result: {:-19} vs {:.16e}",
                    format!(
                        "{:-19}",
                        utils::format_uncertainty(res.result[0], res.error[0])
                    )
                    .red()
                    .bold(),
                    target.re
                );
                return false;
            }
            (*settings).integrator.integrated_phase = IntegratedPhase::Imag;
            let res = havana_integrate(&settings, user_data_generator, Some(target));
            if !approx_eq(res.result[1], target.im, applied_tolerance)
                || !validate_error(res.error[1], target.re - res.result[1])
            {
                println!(
                    "Incorrect imag part of result: {:-19} vs {:.16e}",
                    format!(
                        "{:-19}",
                        utils::format_uncertainty(res.result[1], res.error[1])
                    )
                    .red()
                    .bold(),
                    target.im
                );
                return false;
            }
        }
        IntegratedPhase::Real => {
            (*settings).integrator.integrated_phase = IntegratedPhase::Real;
            let res = havana_integrate(&settings, user_data_generator, Some(target));
            if !approx_eq(res.result[0], target.re, applied_tolerance)
                || !validate_error(res.error[0], target.im - res.result[0])
            {
                println!(
                    "Incorrect real part of result: {:-19} vs {:.16e}",
                    format!(
                        "{:-19}",
                        utils::format_uncertainty(res.result[0], res.error[0])
                    )
                    .red()
                    .bold(),
                    target.re
                );
                return false;
            }
        }
        IntegratedPhase::Imag => {
            (*settings).integrator.integrated_phase = IntegratedPhase::Imag;
            let res = havana_integrate(&settings, user_data_generator, Some(target));
            if !approx_eq(res.result[1], target.im, applied_tolerance)
                || !validate_error(res.error[1], target.im - res.result[1])
            {
                println!(
                    "Incorrect imag part of result: {:-19} vs {:.16e}",
                    format!(
                        "{:-19}",
                        utils::format_uncertainty(res.result[1], res.error[1])
                    )
                    .red()
                    .bold(),
                    target.im
                );
                return false;
            }
        }
    }
    true
}

fn compare_inspect(
    settings: &mut Settings,
    pt: Vec<f64>,
    is_momentum_space: bool,
    target: Complex<f64>,
) -> bool {
    let mut integrand = integrand_factory(&settings);
    let res = inspect(
        &settings,
        &mut integrand,
        pt,
        false,
        is_momentum_space,
        true,
    );
    if !approx_eq(res.re, target.re, INSPECT_TOLERANCE)
        || !approx_eq(res.im, target.im, INSPECT_TOLERANCE)
    {
        println!(
            "Incorrect result from inspect: {}\n                            vs {}",
            format!("{:+16e} + i {:+16e}", res.re, res.im).red().bold(),
            format!("{:.16e} + i {:+16e}", target.re, target.im)
                .red()
                .bold()
        );
        return false;
    }
    return true;
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
    use super::*;

    #[test]
    fn unit_volume_11_momenta_hyperspherical_flat() {
        let mut settings = load_default_settings();
        let mut itg = get_unit_volume_integrand();
        settings.integrator.n_start = 5 * BASE_N_START_SAMPLE;
        settings.integrator.n_max = 10 * BASE_N_START_SAMPLE;
        settings.integrator.n_increase = 0;
        settings.integrator.n_increase = 0;
        settings.kinematics.e_cm = 1.;

        itg.n_3d_momenta = 11;
        settings.parameterization.mode = ParameterizationMode::HyperSphericalFlat;
        settings.hard_coded_integrand = IntegrandSettings::UnitVolume(itg.clone());
        assert!(compare_integration(
            &mut settings,
            IntegratedPhase::Real,
            Complex::new(1.0, 0.0),
            None
        ));
    }

    #[test]
    fn unit_volume_3_momenta_hyperspherical() {
        let mut settings = load_default_settings();
        let mut itg = get_unit_volume_integrand();
        settings.integrator.n_start = 5 * BASE_N_START_SAMPLE;
        settings.integrator.n_max = 10 * BASE_N_START_SAMPLE;
        settings.integrator.n_increase = 0;
        settings.integrator.n_increase = 0;
        settings.kinematics.e_cm = 1.;

        itg.n_3d_momenta = 3;
        settings.parameterization.mode = ParameterizationMode::HyperSpherical;
        settings.hard_coded_integrand = IntegrandSettings::UnitVolume(itg.clone());
        assert!(compare_integration(
            &mut settings,
            IntegratedPhase::Real,
            Complex::new(1.0, 0.0),
            None
        ));
    }

    #[test]
    fn unit_volume_3_momenta_spherical() {
        let mut settings = load_default_settings();
        let mut itg = get_unit_volume_integrand();
        settings.integrator.n_start = 5 * BASE_N_START_SAMPLE;
        settings.integrator.n_max = 10 * BASE_N_START_SAMPLE;
        settings.integrator.n_increase = 0;
        settings.integrator.n_increase = 0;
        settings.kinematics.e_cm = 1.;

        itg.n_3d_momenta = 3;
        settings.parameterization.mode = ParameterizationMode::Spherical;
        settings.hard_coded_integrand = IntegrandSettings::UnitVolume(itg.clone());
        assert!(compare_integration(
            &mut settings,
            IntegratedPhase::Real,
            Complex::new(1.0, 0.0),
            None
        ));
    }

    #[test]
    fn poly_left_right_exponential_h_function() {
        let mut settings = load_default_settings();
        let mut itg = get_h_function_test_integrand();
        settings.integrator.n_start = 5 * BASE_N_START_SAMPLE;
        settings.integrator.n_max = 10 * BASE_N_START_SAMPLE;
        settings.integrator.n_increase = 0;
        settings.integrator.n_increase = 0;
        settings.kinematics.e_cm = 1.;

        itg.h_function = HFunctionSettings {
            function: HFunction::PolyLeftRightExponential,
            sigma: 0.01,
            power: Some(12),
            enabled_dampening: true,
        };
        settings.hard_coded_integrand = IntegrandSettings::HFunctionTest(itg.clone());
        assert!(compare_integration(
            &mut settings,
            IntegratedPhase::Real,
            Complex::new(1.0, 0.0),
            None
        ));
    }
}

#[cfg(test)]
mod tests_inspect {
    use super::*;

    // Amazingly enough, a simple failing test induces a segfault on MacOS... :/
    // #[test]
    // fn this_test_will_not_pass() {
    //     assert!(false);
    // }

    #[test]
    fn inspect_unit_volume() {
        let mut settings = load_default_settings();
        let mut itg = get_unit_volume_integrand();
        itg.n_3d_momenta = 6;
        settings.kinematics.e_cm = 1.;
        settings.parameterization.mode = ParameterizationMode::Spherical;
        settings.hard_coded_integrand = IntegrandSettings::UnitVolume(itg.clone());
        assert!(compare_inspect(
            &mut settings,
            vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            true,
            Complex::new(1.3793965770302298e3, 0.0)
        ));

        itg.n_3d_momenta = 9;
        settings.kinematics.e_cm = 100.;
        settings.parameterization.mode = ParameterizationMode::HyperSpherical;
        settings.hard_coded_integrand = IntegrandSettings::UnitVolume(itg.clone());
        assert!(compare_inspect(
            &mut settings,
            vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            true,
            Complex::new(4.792927924134406e-45, 0.0)
        ));
    }

    #[test]
    fn inspect_h_function_test() {
        let mut settings = load_default_settings();
        let mut itg = get_h_function_test_integrand();
        settings.kinematics.e_cm = 1.;

        itg.h_function = HFunctionSettings {
            function: HFunction::PolyLeftRightExponential,
            sigma: 0.01,
            power: Some(12),
            enabled_dampening: true,
        };
        settings.hard_coded_integrand = IntegrandSettings::HFunctionTest(itg.clone());
        assert!(compare_inspect(
            &mut settings,
            vec![0.2188450233532342,],
            false,
            Complex::new(1.4016882047579115e-34, 0.0)
        ));

        itg.h_function = HFunctionSettings {
            function: HFunction::PolyLeftRightExponential,
            sigma: 0.3,
            power: Some(9),
            enabled_dampening: false,
        };
        settings.hard_coded_integrand = IntegrandSettings::HFunctionTest(itg.clone());
        assert!(compare_inspect(
            &mut settings,
            vec![0.2188450233532342,],
            false,
            Complex::new(3.112977432926161e-4, 0.0)
        ));
    }
}
