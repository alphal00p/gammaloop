use enum_dispatch::enum_dispatch;

use crate::evaluation_result::EvaluationResult;
use crate::integrands::{HasIntegrand, Integrand};
use crate::integrate::UserData;
use crate::utils::F;
use symbolica::numerical_integration::{Grid, Sample};
pub mod amplitude_integrand;
pub mod cross_section_integrand;
use crate::observables::EventManager;
use crate::{IntegratorSettings, Settings};

#[derive(Clone)]
#[enum_dispatch(HasIntegrand)]
pub enum NewIntegrand {
    Amplitude(amplitude_integrand::AmplitudeIntegrand),
    CrossSection(cross_section_integrand::CrossSectionIntegrand),
}

impl NewIntegrand {
    pub fn get_settings(&self) -> &Settings {
        match self {
            NewIntegrand::Amplitude(integrand) => &integrand.settings,
            NewIntegrand::CrossSection(integrand) => &integrand.settings,
        }
    }

    /// Used to create the use_data_generator closure for havana_integrate
    pub fn user_data_generator(&self, num_cores: usize, _settings: &Settings) -> UserData {
        UserData {
            integrand: vec![Integrand::NewIntegrand(self.clone()); num_cores],
        }
    }

    pub fn get_mut_settings(&mut self) -> &mut Settings {
        match self {
            NewIntegrand::Amplitude(integrand) => &mut integrand.settings,
            NewIntegrand::CrossSection(integrand) => &mut integrand.settings,
        }
    }
}
