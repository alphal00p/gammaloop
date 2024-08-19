use std::sync::{LazyLock, Mutex, MutexGuard};

use gammaloop_integrand::GammaLoopSample;
use serde::{Deserialize, Serialize};
use symbolica::numerical_integration::Sample;

use crate::{gammaloop_integrand, utils::F};

pub static DEBUG_INFO: DebugInfo = DebugInfo::new();

pub struct DebugInfo {
    data: LazyLock<Mutex<DebugInfoData>>,
}

impl DebugInfo {
    const fn new() -> Self {
        Self {
            data: LazyLock::new(|| Mutex::new(DebugInfoData::new())),
        }
    }

    #[cold]
    pub fn get(&self) -> MutexGuard<'_, DebugInfoData> {
        self.data.lock().unwrap()
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct DebugInfoData {
    pub havana_sample: Option<Sample<F<f64>>>,
    pub gammaloop_sample: Option<GammaLoopSample<f64>>,
}

impl DebugInfoData {
    fn new() -> Self {
        Self {
            havana_sample: None,
            gammaloop_sample: None,
        }
    }
}
