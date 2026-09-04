use std::ops::Neg;

use color_eyre::Result;
use eyre::eyre;
#[cfg(test)]
use symbolica::atom::Atom;

use crate::{
    cff::{expression::OrientationID, surface::LinearEnergyExpr},
    graph::{LoopMomentumBasis, cuts::CutSet},
    uv::{
        Integrands,
        approx::{OrientationProjection, direct_3d::Direct3dCts, projected_4d::Projected4dCts},
    },
};

mod integrated_localizer;
mod projection_branches;
mod residue_localizer;

#[cfg(test)]
mod tests;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct FrozenActiveCt {
    pub active: OrientationIntegrands,
    pub frozen_integrands: Integrands,
    /// The coordinate frame in which this sector's still-active Taylor
    /// coefficient was formed. Direct complete-CFF sectors have no such
    /// independent 4D frame.
    pub active_lmb: Option<LoopMomentumBasis>,
}

/// Residue integrands grouped by the selector and exact energy map that own
/// every numerator factor in the term. A missing source map uses the selected
/// production map; a present source map remains authoritative while production
/// IDs only partition exact residue-map-key hosts.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(super) struct OrientationIntegrandBranch {
    pub(super) selector_id: OrientationID,
    pub(super) source_edge_energy_map: Option<Vec<LinearEnergyExpr>>,
    pub(super) integrands: Integrands,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct OrientationIntegrands(pub(super) Vec<OrientationIntegrandBranch>);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SourceSelectorHosting {
    /// Direct 3D branches inherit one compatible complete production key.
    PhysicalPrefix,
    /// Projected local-4D sources are summed independently. Prefer a compatible
    /// production key as stable mapping metadata, but fall back to any
    /// permitted deterministic host when a generalized source has no complete
    /// physical extension. The source-local energy map owns evaluation.
    IndependentSum,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Localizer<'a> {
    pub(super) cutset: &'a CutSet,
    pub(super) orientation: OrientationProjection<'a>,
    pub(super) source_selector_hosting: SourceSelectorHosting,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Local3DCts {
    /// The typed direct lane. Its complete generalized residue keys remain
    /// opaque while the Taylor forest is built.
    Direct(Direct3dCts),
    Projected4d(Projected4dCts),
}

impl Neg for Local3DCts {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            Self::Direct(direct) => Self::Direct(-direct),
            Self::Projected4d(sectors) => Self::Projected4d(-sectors),
        }
    }
}

impl Local3DCts {
    pub(crate) fn direct(&self) -> Result<&Direct3dCts> {
        match self {
            Self::Direct(direct) => Ok(direct),
            _ => Err(eyre!(
                "projected local-4D coefficients are not direct local-3D counterterms"
            )),
        }
    }

    #[cfg(test)]
    pub(crate) fn map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        match self {
            Self::Direct(direct) => Ok(Self::Direct(direct.map(&mut f)?)),
            Self::Projected4d(sectors) => Ok(Self::Projected4d(sectors.map(&mut f)?)),
        }
    }
}
