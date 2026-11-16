use std::collections::{BTreeMap, BTreeSet};

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubSetLike};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;
use three_dimensional_reps::{MediumMode, ThermalDistributionCall, ThermalDistributionFactor};

use crate::graph::LoopMomentumBasis;
use crate::integrands::process::param_builder::ParamBuilderGraph;
use crate::utils::{GS, W_};
use crate::uv::approx::UVCtx;
use crate::uv::uv_graph::UVE;

use super::{
    Direct3dApproximation, DirectResidueBranches, LOCAL_3D_SOFT_SCOPE, Local3DDeformation,
};

/// Decode a component owner written as the `S_<label>` symbol of its support.
pub(super) fn component_owner(owner: AtomView<'_>, n_hedges: usize) -> Option<SuBitGraph> {
    let AtomView::Var(owner) = owner else {
        return None;
    };
    owner
        .get_symbol()
        .get_stripped_name()
        .strip_prefix("S_")
        .and_then(|label| SuBitGraph::from_base62(label, n_hedges))
}

/// Mass-scale tag of the completed thermal soft branch of `support`.
pub(super) fn soft_branch_scope(support: &SuBitGraph) -> Atom {
    function!(*LOCAL_3D_SOFT_SCOPE, support.symbol())
}

/// Relation of a completed thermal soft branch to one local operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum SoftBranchRole {
    /// Owned by a strictly nested soft operation: a frozen vertex whose
    /// scales are co-scaled with the operation's loop momenta.
    Frozen,
    /// This component's own branch, acted on by its own U.
    Own,
    /// Outside the rescaled support, so left unchanged.
    Disjoint,
}

/// Completed thermal soft branches as seen by one local operation.
///
/// A branch carries its owner in two places: its masses are multiplied by
/// [`soft_branch_scope`], and every on-shell energy it expanded uses the
/// owner's `S_<label>` symbol in place of an integer route owner.
#[derive(Clone, Copy)]
pub(super) struct SoftBranches<'a> {
    support: &'a SuBitGraph,
    n_hedges: usize,
    acts_on_own_branch: bool,
}

impl<'a> SoftBranches<'a> {
    pub(super) fn new(
        support: &'a SuBitGraph,
        n_hedges: usize,
        deformation: Local3DDeformation,
    ) -> Self {
        Self {
            support,
            n_hedges,
            acts_on_own_branch: deformation == Local3DDeformation::OrdinaryOnSoftBranch,
        }
    }

    pub(super) fn role(&self, owner: &SuBitGraph) -> Result<SoftBranchRole> {
        if owner == self.support {
            // Only the component's own U revisits the support of its soft branch.
            return if self.acts_on_own_branch {
                Ok(SoftBranchRole::Own)
            } else {
                Err(eyre!(
                    "the completed soft branch of {} reached another operation on the same support",
                    owner.string_label()
                ))
            };
        }
        if self.support.includes(owner) {
            Ok(SoftBranchRole::Frozen)
        } else if self.support.intersects(owner) {
            Err(eyre!(
                "the completed soft branch of {} partially overlaps the operation on {}",
                owner.string_label(),
                self.support.string_label()
            ))
        } else {
            Ok(SoftBranchRole::Disjoint)
        }
    }

    /// The role of a soft-branch mass scope or energy owner, if `view` is one.
    pub(super) fn tag_role(&self, view: AtomView<'_>) -> Result<Option<SoftBranchRole>> {
        let AtomView::Fun(call) = view else {
            return Ok(None);
        };
        let owner = if call.get_symbol() == *LOCAL_3D_SOFT_SCOPE {
            (call.get_nargs() == 1)
                .then(|| component_owner(call.get(0), self.n_hedges))
                .flatten()
                .ok_or_else(|| eyre!("a thermal soft-branch scope has no valid owner: {view}"))?
        } else if call.get_symbol() == GS.on_shell_energy && call.get_nargs() == 2 {
            // Integer route owners belong to ordinary operations.
            let Some(owner) = component_owner(call.get(0), self.n_hedges) else {
                return Ok(None);
            };
            owner
        } else {
            return Ok(None);
        };
        self.role(&owner).map(Some)
    }

    /// Whether `atom` depends on a frozen soft branch.
    pub(super) fn contains_frozen(&self, atom: AtomView<'_>) -> Result<bool> {
        let mut frozen = false;
        let mut error = None;
        atom.visitor(&mut |view| {
            if frozen || error.is_some() {
                return false;
            }
            match self.tag_role(view) {
                Ok(Some(SoftBranchRole::Frozen)) => {
                    frozen = true;
                    false
                }
                Ok(Some(_)) | Ok(None) => true,
                Err(tag_error) => {
                    error = Some(tag_error);
                    false
                }
            }
        });
        error.map_or(Ok(frozen), Err)
    }
}

impl Direct3dApproximation<'_> {
    /// The vacuum limit an ordinary operation takes of the occupations on its
    /// support. Occupations retained by a frozen soft branch are kept: their
    /// temperature is co-scaled by [`Self::compensate_thermal_distributions`].
    pub(super) fn project_unfrozen_occupations(
        ctx: &UVCtx<'_>,
        integrand: &Atom,
        support: &SuBitGraph,
        deformation: Local3DDeformation,
    ) -> Result<Atom> {
        let soft_branches = SoftBranches::new(support, ctx.graph.n_hedges(), deformation);
        let active_edges = ctx
            .graph
            .iter_edges_of(support)
            .map(|(_, edge, _)| edge)
            .collect::<BTreeSet<_>>();
        ctx.graph.make_selected_thermal_distributions_explicit(
            integrand,
            MediumMode::Vacuum,
            |call| {
                Ok(active_edges.contains(&call.edge)
                    && !soft_branches.contains_frozen(call.energy.as_view())?)
            },
        )
    }

    /// Hand a completed thermal soft branch to enclosing operations as a frozen
    /// vertex. Its mass scope becomes the branch scope, and every on-shell energy
    /// this operation expanded becomes owned by the support's symbol. Branches
    /// that were already frozen inside it keep their own owners.
    pub(super) fn freeze_soft_branch(
        branch: &DirectResidueBranches,
        support: &SuBitGraph,
        lmb: &LoopMomentumBasis,
        mass_scope: &Atom,
    ) -> Result<DirectResidueBranches> {
        branch.map_expressions(|atom| Ok(Self::freeze_soft_atom(atom, support, lmb, mass_scope)))
    }

    /// [`Self::freeze_soft_branch`] for one atom.
    pub(super) fn freeze_soft_atom(
        atom: &Atom,
        support: &SuBitGraph,
        lmb: &LoopMomentumBasis,
        mass_scope: &Atom,
    ) -> Atom {
        let owner = Atom::var(support.symbol());
        let mut frozen = atom
            .replace(mass_scope.clone())
            .with(soft_branch_scope(support));
        for edge in &lmb.loop_edges {
            frozen = frozen
                .replace(function!(
                    GS.on_shell_energy,
                    usize::from(*edge) as i64,
                    W_.prop_
                ))
                .with(function!(GS.on_shell_energy, owner.clone(), W_.prop_));
        }
        frozen
    }

    /// Express the soft chart at fixed physical temperature and chemical potential.
    /// For ordinary energy derivatives, N_d(h E; h T, h mu) = h^-d N_d(E; T, mu).
    /// The routed energy is h E(k + q/h, m), so dividing its argument by h
    /// exposes the external-momentum Taylor series. The factor h^-d retains
    /// the derivative's engineering dimension in the existing Laurent counting.
    ///
    /// The soft deformation compensates every occupation of its active
    /// component. Any deformation also compensates the occupations retained by
    /// a frozen soft branch, keeping their owners, since that branch is a vertex
    /// whose temperature and chemical potential scale with its loop momenta.
    pub(super) fn compensate_thermal_distributions(
        ctx: &UVCtx<'_>,
        active_subgraph: &SuBitGraph,
        integrand: &Atom,
        deformation: Local3DDeformation,
        soft_branches: SoftBranches<'_>,
    ) -> Result<Atom> {
        // N is the same function for equal statistics and chemical potential.
        // Choose its representative within this active component so coincident
        // energies cancel symbolically. Ownership of the branch lives in the
        // energy's owner, not in this representative; disjoint components and
        // their occupations remain unchanged.
        let mut representatives = BTreeMap::new();
        let owners = if deformation == Local3DDeformation::Soft {
            ctx.graph
                .iter_edges_of(active_subgraph)
                .filter_map(|(pair, edge, _)| {
                    if !pair.is_paired() {
                        return None;
                    }
                    let identity = (
                        ctx.graph[edge].is_fermion(),
                        ctx.graph[edge].chemical_potential_atom(),
                    );
                    Some((edge, *representatives.entry(identity).or_insert(edge)))
                })
                .collect::<BTreeMap<_, _>>()
        } else {
            BTreeMap::new()
        };
        let scale = Atom::var(GS.rescale);
        let mut error = None;
        let compensated = integrand.replace_map(|view, _, output| {
            if error.is_some() {
                return;
            }
            let AtomView::Fun(call) = view else {
                return;
            };
            if call.get_symbol() != GS.thermal_distribution {
                return;
            }
            let distribution = match ThermalDistributionCall::try_from(view) {
                Ok(distribution) => distribution,
                Err(parse_error) => {
                    error = Some(parse_error.into());
                    return;
                }
            };
            if !distribution.temperature_flag.is_one() {
                return;
            }
            let owner = match soft_branches.contains_frozen(distribution.energy.as_view()) {
                Ok(true) => distribution.edge,
                Ok(false) => match owners.get(&distribution.edge) {
                    Some(&owner) => owner,
                    None => return,
                },
                Err(tag_error) => {
                    error = Some(tag_error);
                    return;
                }
            };
            let vacuum = ThermalDistributionFactor::vacuum_atom(
                distribution.derivative_order,
                distribution.thermal_sign,
            );
            let derivative_scale = scale.pow(-Atom::num(distribution.derivative_order));
            // N_0(s) = (1+s)/2 + N_0(-1); positive derivatives do
            // not depend on s. One occupation contour makes cyclic-chain
            // derivatives and Taylor derivatives identical.
            **output = vacuum
                + derivative_scale
                    * ThermalDistributionCall::atom_at_energy(
                        owner.0,
                        distribution.derivative_order,
                        distribution.temperature_flag,
                        -1,
                        distribution.orientation_sign,
                        distribution.energy / &scale,
                    );
        });
        error.map_or(Ok(compensated), Err)
    }
}
