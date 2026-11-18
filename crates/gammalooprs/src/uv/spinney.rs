use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

use color_eyre::Result;
use eyre::WrapErr;
use gammaloop_tracing_filter::LogMessage;
use linnet::half_edge::{
    HedgeGraph,
    subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike},
};
use tracing::debug;

use crate::debug_tags;
use crate::uv::thermal_zero_mode::{ThermalZeroModeRegion, ZeroModeHostKind};
use crate::{
    graph::{LMBext, LoopMomentumBasis, cuts::CutSet},
    uv::{ApproximationType, UltravioletGraph},
};

#[derive(Clone, Debug, Eq)]
pub struct Spinney {
    pub subgraph: InternalSubGraph,
    components: Vec<SuBitGraph>,
    pub dod: i32,
    pub lmb: LoopMomentumBasis,
    /// Disconnected unions use `MUV` as a neutral tag; their component factors
    /// carry the renormalization schemes that determine the approximation.
    /// A union that hosts a thermal zero-mode region is the exception: its
    /// copy carries `IR` and the region.
    pub renormalization_scheme: ApproximationType,
    max_comp_loop_count: usize,
    /// The finite-temperature zero-mode region this spinney hosts, if any.
    thermal_zero_mode: Option<ThermalZeroModeRegion>,
}

/// The node a spinney becomes in the hedge-poset wood. A union host shares
/// its subgraph with the plain factorized join of its components, and sits
/// above it. The wood keys its nodes on the subgraph and this role, and
/// treats a union host as one collective operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum WoodRole {
    Plain,
    UnionHost,
}

impl LogMessage for Spinney {
    fn log_display(&self) -> String {
        let region = self
            .thermal_zero_mode
            .as_ref()
            .map(|region| format!(", thermal_zero_mode: {region}"))
            .unwrap_or_default();
        format!(
            "Spinney {{ dod: {}, components: {}, max_comp_loop_count: {}, scheme: {}, subgraph: {}{region} }}",
            self.dod,
            self.components.len(),
            self.max_comp_loop_count,
            self.renormalization_scheme,
            self.subgraph.string_label(),
        )
    }
}

impl Spinney {
    pub fn compatible_with(&self, cut: &CutSet) -> bool {
        !self.subgraph.filter.intersects(&cut.union)
    }

    pub fn empty<E, V, H, G: AsRef<HedgeGraph<E, V, H>> + LMBext + ?Sized>(g: &G) -> Self {
        Self {
            subgraph: InternalSubGraph::empty(g.as_ref().n_hedges()),
            components: vec![],
            dod: 0,
            lmb: g.empty_lmb(),
            renormalization_scheme: ApproximationType::MUV,
            max_comp_loop_count: 0,
            thermal_zero_mode: None,
        }
    }

    pub fn new<E, V, H, G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>> + ?Sized>(
        subgraph: InternalSubGraph,
        g: &G,
        lmb: &LoopMomentumBasis,
    ) -> Result<Option<Self>> {
        let dod = g.compute_dod(&subgraph);
        Self::with_scheme(subgraph, g, lmb, ApproximationType::MUV, dod)
    }

    pub fn with_scheme<E, V, H, G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>> + ?Sized>(
        subgraph: InternalSubGraph,
        g: &G,
        lmb: &LoopMomentumBasis,
        renormalization_scheme: ApproximationType,
        dod: i32,
    ) -> Result<Option<Self>> {
        let components = g.as_ref().connected_components(&subgraph);
        let max_comp_loop_count = components
            .iter()
            .map(|component| g.as_ref().cyclotomatic_number(component))
            .max()
            .unwrap_or(0);

        // Divergence is a property of each connected factor of a spinney.
        // A disconnected union may have a negative aggregate DOD even though
        // every factor was selected as divergent and must still be retained.
        if dod < 0 && components.len() <= 1 {
            return Ok(None);
        }

        let lmb = g
            .try_compatible_sub_lmb(&subgraph, g.dummy_less_full_crown(&subgraph), lmb)
            .wrap_err_with(|| {
                format!(
                    "failed to select a compatible loop-momentum route for spinney {} (scheme {renormalization_scheme}, DOD {dod})",
                    subgraph.string_label()
                )
            })?;

        debug!(
            dod = %dod,
            graph = %g.dot_lmb_of(&subgraph,&lmb),
            renormalization_scheme = %renormalization_scheme,
            string_label =%subgraph.string_label(),
            "Constucted spinney with scheme: {:?}",
            renormalization_scheme
        );

        Ok(Some(Self {
            components,
            dod,
            lmb,
            subgraph,
            renormalization_scheme,
            max_comp_loop_count,
            thermal_zero_mode: None,
        }))
    }

    /// This spinney as the host of `region`. A host is subtracted by the
    /// zero-mode operation, so its scheme is `IR`; a union host is a copy of
    /// the plain union, which keeps its neutral scheme.
    pub(crate) fn with_thermal_zero_mode(mut self, region: ThermalZeroModeRegion) -> Self {
        self.renormalization_scheme = ApproximationType::IR;
        self.thermal_zero_mode = Some(region);
        debug_tags!(#uv, #generation, #zero_mode;
            log.spinney = self,
            "Assigned a thermal zero-mode region"
        );
        self
    }

    /// This spinney as a plain `MUV` U. A zero-mode host is subtracted by its
    /// zero-mode operation and keeps `IR`.
    ///
    /// # Panics
    ///
    /// If this spinney hosts a zero-mode region: its callers demote only
    /// spinneys that host none.
    pub(crate) fn demote_to_muv(&mut self) {
        assert!(
            self.thermal_zero_mode.is_none(),
            "the zero-mode host {} keeps its IR scheme",
            self.subgraph.string_label()
        );
        self.renormalization_scheme = ApproximationType::MUV;
    }

    /// The finite-temperature zero-mode region this spinney hosts, if any.
    pub(crate) fn thermal_zero_mode(&self) -> Option<&ThermalZeroModeRegion> {
        self.thermal_zero_mode.as_ref()
    }

    pub fn filter(&self) -> &SuBitGraph {
        &self.subgraph.filter
    }

    pub fn max_comp_loop_count(&self) -> usize {
        self.max_comp_loop_count
    }

    pub fn n_components(&self) -> usize {
        self.components.len()
    }

    /// The connected components of the subgraph.
    pub(crate) fn components(&self) -> &[SuBitGraph] {
        &self.components
    }

    /// The node this spinney becomes in the hedge-poset wood.
    pub(crate) fn wood_role(&self) -> WoodRole {
        match self
            .thermal_zero_mode
            .as_ref()
            .map(ThermalZeroModeRegion::kind)
        {
            Some(ZeroModeHostKind::Union) => WoodRole::UnionHost,
            Some(ZeroModeHostKind::Connected) | None => WoodRole::Plain,
        }
    }
}

impl Hash for Spinney {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.subgraph.hash(state);
        self.wood_role().hash(state);
    }
}

/// Spinneys are equal when they are the same wood node: the same subgraph
/// in the same role.
impl PartialEq for Spinney {
    fn eq(&self, other: &Self) -> bool {
        self.subgraph == other.subgraph && self.wood_role() == other.wood_role()
    }
}

/// Inclusion of subgraphs; for one subgraph, the plain node precedes its
/// union host.
impl PartialOrd for Spinney {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.subgraph == other.subgraph {
            Some(self.wood_role().cmp(&other.wood_role()))
        } else if self.subgraph.includes(&other.subgraph) {
            Some(Ordering::Greater)
        } else if other.subgraph.includes(&self.subgraph) {
            Some(Ordering::Less)
        } else {
            None
        }
    }
}
