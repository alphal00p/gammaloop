//! `oneloop`: one-loop IBP reduction

pub mod error;
pub mod family;
pub mod masters;
pub mod reduce;
pub mod solver;
pub mod symbols;

pub use error::OneLoopError;
pub use family::{Integral, IntegralFamily, Isp, Kinematics, Propagator};
pub use masters::{MasterBasis, MasterIntegral, OneLoopMasters};
pub use reduce::{Reduction, reduce};
pub use solver::RationalSolver;

pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
pub(crate) fn ensure_symbolica_license() {
    use std::sync::Once;
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        gammalooprs::initialisation::initialise()
            .expect("gammaloop initialisation (activates the Symbolica license)");
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_and_reports_version() {
        assert_eq!(version(), "0.1.0");
    }
}
