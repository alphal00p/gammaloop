//! `oneloop`: one-loop IBP reduction

pub mod amplitude;
pub mod error;
pub mod family;
pub mod masters;
pub mod reduce;
pub mod symbols;

pub use amplitude::amplitude;
pub use error::OneLoopError;
pub use family::{Integral, IntegralFamily, Isp, Kinematics, Propagator};
pub use masters::{MasterBasis, MasterIntegral, OneLoopMasters};
pub use reduce::{Reduction, reduce};

#[cfg(test)]
pub(crate) fn ensure_symbolica_license() {
    use std::sync::Once;
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        gammalooprs::initialisation::initialise()
            .expect("gammaloop initialisation (activates the Symbolica license)");
    });
}
