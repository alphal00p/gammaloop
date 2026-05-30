//! `oneloop`: one-loop IBP reduction

pub mod error;
pub mod family;
pub mod masters;

pub use error::OneLoopError;
pub use family::{Integral, IntegralFamily, Isp, Kinematics, Propagator};
pub use masters::{MasterBasis, MasterIntegral, OneLoopMasters};

pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Initialise gammaloop once per test process. This activates the workspace's
/// Symbolica license (from within `gammalooprs`, where the OEM key validates),
/// so multiple Symbolica-using tests can share a single process. Test-only —
/// the library itself never initialises a license; that is the consuming
/// binary's responsibility.
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
