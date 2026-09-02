use std::{collections::BTreeMap, ops::Neg};

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::subgraph::{SuBitGraph, SubSetOps};
use symbolica::atom::Atom;

use crate::{
    cff::{CutCFFIndex, expression::OrientationID},
    graph::Graph,
    utils::GS,
    uv::{
        Integrands, UVgenerationSettings, UltravioletGraph,
        approx::{
            ForestNodeLike, UVCtx,
            integrated::IntegratedCts,
            local_3d::{FrozenActiveCt, Localizer},
        },
        marker::{UvMarker, UvOperation},
    },
};

use super::{
    branches::{DirectResidueBranches, DirectResidueKey},
    kernel::{Local3DLoopRescaling, apply_taylor},
};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct DirectSector {
    pub(crate) active_subgraph: SuBitGraph,
    pub(crate) active: DirectResidueBranches,
    pub(crate) frozen_integrands: Integrands,
}

impl DirectSector {
    pub(crate) fn combine(&self) -> Result<DirectResidueBranches> {
        self.active.zip_mul_unmapped(&self.frozen_integrands)
    }
}

impl Neg for DirectSector {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            active_subgraph: self.active_subgraph,
            active: -self.active,
            frozen_integrands: self.frozen_integrands,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Direct3dCts {
    /// The empty forest owns the complete stored production CFF.
    Root(DirectResidueBranches),
    /// Direct local-3D forest sectors. Every active factor contains the
    /// complete post-energy-integration CFF on which subsequent Taylor
    /// operators act.
    Sectors(Vec<DirectSector>),
}

impl Neg for Direct3dCts {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            Self::Root(branches) => Self::Root(-branches),
            Self::Sectors(sectors) => {
                Self::Sectors(sectors.into_iter().map(DirectSector::neg).collect())
            }
        }
    }
}

impl Direct3dCts {
    /// Load the exact stored production CFF without passing through reduced-
    /// source representative reconstruction. The production ID binds the
    /// complete loop- and edge-energy map recoverable from the production
    /// expression. The full branch key remains an opaque sparse selector while
    /// Taylor operators act; no selector derived from coarser physical edge
    /// directions is materialized here.
    pub(crate) fn root(graph: &Graph, localizer: Localizer<'_>) -> Result<Self> {
        let production = localizer.orientation.root_expression().ok_or_else(|| {
            eyre!("the direct local-3D root requires its stored production CFF expression")
        })?;
        let cff = graph.cff_from_production_expression(
            production,
            localizer.cutset,
            localizer.orientation.orientation_pattern,
        )?;
        localizer
            .orientation
            .record_energy_degree_bound_report(&cff.energy_degree_bound_report);
        let indices = cff.terms.keys().copied().collect::<Vec<_>>();
        let production_prefactor = Atom::num(cff.production_prefactor_factor());
        let exact_orientations = localizer
            .orientation
            .exact_orientations()
            .expect("a stored production expression has exact residue maps");
        let mut branches: BTreeMap<OrientationID, BTreeMap<CutCFFIndex, Atom>> = BTreeMap::new();

        for (index, term) in cff.terms {
            for stored in term.orientations {
                let id = stored.production_orientation_id.ok_or_else(|| {
                    eyre!("a stored production CFF entry has no production orientation ID")
                })?;
                let original = exact_orientations
                    .get(id)
                    .ok_or_else(|| eyre!("missing stored production orientation {}", id.0))?;
                if stored.orientation.data.orientation != original.data.orientation
                    || stored.orientation.loop_energy_map != original.loop_energy_map
                    || stored.orientation.edge_energy_map != original.edge_energy_map
                {
                    return Err(eyre!(
                        "stored production CFF entry {} changed its exact residue map",
                        id.0
                    ));
                }
                let branch = branches
                    .entry(id)
                    .or_insert_with(|| indices.iter().map(|index| (*index, Atom::Zero)).collect());
                *branch
                    .get_mut(&index)
                    .expect("all direct-root CutCFF indices were initialized") +=
                    stored.expression * &production_prefactor;
            }
        }
        if branches.is_empty() {
            return Err(eyre!(
                "orientation pattern selects no production residue maps"
            ));
        }

        let four_d_denominators = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );
        Ok(Self::Root(DirectResidueBranches::from_keyed(
            branches.into_iter().map(|(id, integrands)| {
                (
                    DirectResidueKey::production(id),
                    Integrands::from_iter(integrands).map(|atom| atom * &four_d_denominators),
                )
            }),
        )?))
    }

    pub(crate) fn from_sectors(sectors: Vec<DirectSector>) -> Result<Self> {
        if sectors.is_empty() {
            return Err(eyre!("active UV counterterm sectors cannot be empty"));
        }
        Ok(Self::Sectors(sectors))
    }

    pub(crate) fn sectors(&self) -> Result<&[DirectSector]> {
        match self {
            Self::Sectors(sectors) => Ok(sectors),
            Self::Root(_) => Err(eyre!("the direct root has not entered a UV sector")),
        }
    }

    /// Materialize only a direct post-energy-integration CFF value. Projected
    /// local-4D coefficients are incomplete until final assembly supplies the
    /// untouched outer CFF and therefore cannot cross this boundary.
    pub(crate) fn branches(&self) -> Result<DirectResidueBranches> {
        match self {
            Self::Root(branches) => Ok(branches.clone()),
            Self::Sectors(sectors) => {
                let mut sectors = sectors.iter();
                let mut branches = sectors
                    .next()
                    .ok_or_else(|| eyre!("active UV counterterm sectors cannot be empty"))?
                    .combine()?;
                for sector in sectors {
                    branches = branches.zip_add(&sector.combine()?)?;
                }
                Ok(branches)
            }
        }
    }

    pub(crate) fn map(&self, mut map: impl FnMut(&Atom) -> Result<Atom>) -> Result<Self> {
        match self {
            Self::Root(branches) => Ok(Self::Root(branches.fallible_map(|_, atom| map(atom))?)),
            Self::Sectors(sectors) => Ok(Self::Sectors(
                sectors
                    .iter()
                    .map(|sector| {
                        Ok(DirectSector {
                            active_subgraph: sector.active_subgraph.clone(),
                            active: sector.active.fallible_map(|_, atom| map(atom))?,
                            frozen_integrands: sector.frozen_integrands.clone(),
                        })
                    })
                    .collect::<Result<_>>()?,
            )),
        }
    }
}

pub(crate) struct Direct3dApproximation<'a> {
    localizer: Localizer<'a>,
    graph: &'a mut Graph,
    settings: &'a UVgenerationSettings,
}

impl<'a> Direct3dApproximation<'a> {
    pub(crate) fn new(
        localizer: Localizer<'a>,
        graph: &'a mut Graph,
        settings: &'a UVgenerationSettings,
    ) -> Self {
        Self {
            localizer,
            graph,
            settings,
        }
    }

    /// Convert the projection helper's transient representation into the
    /// direct lane's explicit full-key type. A source map remains attached to
    /// its branch and therefore governs every later factor in that branch.
    fn localize_integrated<S: ForestNodeLike>(
        &mut self,
        expr: &Atom,
        integrated_node: &S,
    ) -> Result<(DirectResidueBranches, Integrands)> {
        let FrozenActiveCt {
            active,
            frozen_integrands,
            active_lmb,
        } = self.localizer.localize(expr, self.graph, integrated_node)?;
        if active_lmb.is_some() {
            return Err(eyre!(
                "a direct integrated localization unexpectedly retained a projected-4D Taylor LMB"
            ));
        }
        let branches = DirectResidueBranches::from_transient(&active)?;
        Ok((branches, frozen_integrands))
    }

    pub(crate) fn run<S: ForestNodeLike, M: ForestNodeLike>(
        mut self,
        local: &Direct3dCts,
        integrated: &IntegratedCts,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Direct3dCts> {
        // Direct local-3D subtraction acts on the complete CFF after the energy
        // integrations. The explicit-sum option changes only whether each
        // branch's complete residue-map key is materialized as a selector; it
        // must not select a different UV construction.
        let integrated_atom = integrated.physical_finite_counterterm_atom();
        let integrated_t = (!integrated_atom.is_zero())
            .then(|| self.localize_integrated(&integrated_atom, given))
            .transpose()?;
        let ctx = UVCtx::new(self.graph, self.settings);
        let marker = UvMarker::new(ctx.settings);
        let reduced_subgraph = current.reduced_subgraph(given);
        let sectors = match local {
            Direct3dCts::Root(branches) => vec![DirectSector {
                active_subgraph: self.graph.empty_subgraph(),
                active: branches.clone(),
                frozen_integrands: branches.identity_integrands(),
            }],
            Direct3dCts::Sectors(sectors) => sectors.clone(),
        };
        let mut next = Vec::with_capacity(sectors.len() + usize::from(integrated_t.is_some()));

        for sector in sectors {
            let active_subgraph = sector.active_subgraph.union(&reduced_subgraph);
            // A normalized integrated prefix is inert under every later Taylor
            // operation. Transform only its complete active CFF and retain the
            // frozen localizing kernel outside the series.
            let active = -sector.active.fallible_map(|key, atom| {
                apply_taylor(
                    Local3DLoopRescaling::FullSubgraph,
                    &ctx,
                    self.localizer.orientation,
                    current,
                    given,
                    Some(active_subgraph.clone()),
                    key,
                    atom,
                )
            })?;
            next.push(DirectSector {
                active_subgraph,
                active,
                frozen_integrands: sector.frozen_integrands,
            });
        }

        if let Some((active, frozen_integrands)) = integrated_t {
            let active = -active.fallible_map(|key, atom| {
                apply_taylor(
                    Local3DLoopRescaling::ReducedSubgraph,
                    &ctx,
                    self.localizer.orientation,
                    current,
                    given,
                    Some(reduced_subgraph.clone()),
                    key,
                    atom,
                )
            })?;
            next.push(DirectSector {
                active_subgraph: reduced_subgraph,
                active,
                frozen_integrands,
            });
        }

        Direct3dCts::from_sectors(next)?.map(|atom| {
            Ok(marker.apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })
    }

    pub(crate) fn run_local<S: ForestNodeLike, M: ForestNodeLike>(
        self,
        local: &Direct3dCts,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Direct3dCts> {
        let ctx = UVCtx::new(self.graph, self.settings);
        let reduced_subgraph = current.reduced_subgraph(given);
        let sectors = match local {
            Direct3dCts::Root(branches) => vec![DirectSector {
                active_subgraph: self.graph.empty_subgraph(),
                active: branches.clone(),
                frozen_integrands: branches.identity_integrands(),
            }],
            Direct3dCts::Sectors(sectors) => sectors.clone(),
        };
        let next = sectors
            .into_iter()
            .map(|sector| {
                let active_subgraph = sector.active_subgraph.union(&reduced_subgraph);
                // Retain the complete active mask for descendants, but rescale
                // only the component-local part covered by this path.
                let rescaled_subgraph = active_subgraph.intersection(current.subgraph());
                Ok(DirectSector {
                    active_subgraph,
                    active: -sector.active.fallible_map(|key, atom| {
                        apply_taylor(
                            Local3DLoopRescaling::FullSubgraph,
                            &ctx,
                            self.localizer.orientation,
                            current,
                            given,
                            Some(rescaled_subgraph.clone()),
                            key,
                            atom,
                        )
                    })?,
                    frozen_integrands: sector.frozen_integrands,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Direct3dCts::from_sectors(next)?.map(|atom| {
            Ok(UvMarker::new(ctx.settings).apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn run_integrated<S: ForestNodeLike, I: ForestNodeLike, M: ForestNodeLike>(
        mut self,
        integrated: &IntegratedCts,
        integrated_node: &I,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Direct3dCts> {
        let (active, frozen_integrands) = self.localize_integrated(
            &integrated.physical_finite_counterterm_atom(),
            integrated_node,
        )?;
        let ctx = UVCtx::new(self.graph, self.settings);
        let active_subgraph = current.reduced_subgraph(given);
        let active = -active.fallible_map(|key, atom| {
            apply_taylor(
                Local3DLoopRescaling::ReducedSubgraph,
                &ctx,
                self.localizer.orientation,
                current,
                given,
                Some(active_subgraph.clone()),
                key,
                atom,
            )
        })?;
        let sector = DirectSector {
            active_subgraph,
            active,
            frozen_integrands,
        };

        Direct3dCts::from_sectors(vec![sector])?.map(|atom| {
            Ok(UvMarker::new(ctx.settings).apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })
    }
}
