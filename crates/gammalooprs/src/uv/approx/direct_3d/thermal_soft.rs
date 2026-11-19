use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;

use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubSetLike, SubSetOps};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;
use three_dimensional_reps::{MediumMode, ThermalDistributionCall, ThermalDistributionFactor};

use crate::graph::{Graph, LoopMomentumBasis};
use crate::integrands::process::param_builder::{FnMapEntry, ParamBuilderGraph};
use crate::utils::{GS, W_};
use crate::uv::approx::UVCtx;
use crate::uv::uv_graph::UVE;

use super::{
    Direct3dApproximation, DirectResidueBranches, LOCAL_3D_SOFT_SCOPE, Local3DDeformation,
    ZERO_MODE_ORDER, ZERO_MODE_TEMPERATURE, zero_mode_class_representative, zero_mode_temperature,
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

/// Whether `support` contains the frequency class whose representative edge
/// `representative` names.
fn contains_class(graph: &Graph, support: &SuBitGraph, representative: AtomView<'_>) -> bool {
    zero_mode_class_representative(graph, representative)
        .is_some_and(|edge| support.includes(&graph.get_edge_subgraph(edge)))
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
/// owner's `S_<label>` symbol in place of an integer route owner. A completed
/// zero-mode branch also owns the zero-mode temperatures of the classes
/// inside its support and its order tag.
#[derive(Clone, Copy)]
pub(super) struct SoftBranches<'a> {
    graph: &'a Graph,
    support: &'a SuBitGraph,
    deformation: Local3DDeformation,
}

impl<'a> SoftBranches<'a> {
    pub(super) fn new(
        graph: &'a Graph,
        support: &'a SuBitGraph,
        deformation: Local3DDeformation,
    ) -> Self {
        Self {
            graph,
            support,
            deformation,
        }
    }

    pub(super) fn role(&self, owner: &SuBitGraph) -> Result<SoftBranchRole> {
        if owner == self.support {
            // Only the component's own U revisits the support of its soft branch.
            return if self.deformation == Local3DDeformation::OrdinaryOnSoftBranch {
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

    /// The owner of a branch tag, if `view` is one: a soft-branch mass scope,
    /// an energy owned by a branch, or the zero-mode temperature or order tag
    /// of a completed zero-mode branch. The zero-mode temperature of a
    /// frequency class, named by its representative edge, is not yet owned by
    /// a branch.
    fn tag_owner(&self, view: AtomView<'_>) -> Result<Option<SuBitGraph>> {
        let AtomView::Fun(call) = view else {
            return Ok(None);
        };
        let symbol = call.get_symbol();
        let n_hedges = self.graph.n_hedges();
        if symbol == *LOCAL_3D_SOFT_SCOPE || symbol == *ZERO_MODE_ORDER {
            (call.get_nargs() == 1)
                .then(|| component_owner(call.get(0), n_hedges))
                .flatten()
                .map(Some)
                .ok_or_else(|| eyre!("a thermal soft-branch tag has no valid owner: {view}"))
        } else if symbol == *ZERO_MODE_TEMPERATURE {
            let owner = (call.get_nargs() == 1)
                .then(|| call.get(0))
                .ok_or_else(|| eyre!("a zero-mode temperature has no owner: {view}"))?;
            match component_owner(owner, n_hedges) {
                Some(owner) => Ok(Some(owner)),
                // A class temperature, which no branch owns yet.
                None if zero_mode_class_representative(self.graph, owner).is_some() => Ok(None),
                None => Err(eyre!("a zero-mode temperature has no valid owner: {view}")),
            }
        } else if symbol == GS.on_shell_energy && call.get_nargs() == 2 {
            // Integer route owners belong to ordinary operations.
            Ok(component_owner(call.get(0), n_hedges))
        } else {
            Ok(None)
        }
    }

    /// The role of a branch tag, if `view` is one ([`Self::tag_owner`]).
    pub(super) fn tag_role(&self, view: AtomView<'_>) -> Result<Option<SoftBranchRole>> {
        self.tag_owner(view)?
            .map(|owner| self.role(&owner))
            .transpose()
    }

    /// The owners of the completed thermal soft branches that each residue
    /// row of `family` depends on, in row order, through its roots or the
    /// numerator, factor and scalar definitions it calls. The rows of a
    /// family share most of their definitions, so each one is visited once.
    ///
    /// # Errors
    ///
    /// A branch tag without a valid owner.
    pub(super) fn row_owners(
        &self,
        family: &DirectResidueBranches,
    ) -> Result<Vec<BTreeSet<SuBitGraph>>> {
        let mut definition_owners = HashMap::<*const FnMapEntry, BTreeSet<SuBitGraph>>::new();
        family
            .iter_keys()
            .map(|(_, roots)| {
                let mut owners = BTreeSet::new();
                for (_, root) in roots.iter() {
                    self.collect_owners(root.as_view(), &mut owners)?;
                }
                for definition in roots.numerators().iter().chain(roots.scalar_definitions()) {
                    let shared = match definition_owners.entry(Arc::as_ptr(definition)) {
                        Entry::Occupied(entry) => entry.into_mut(),
                        Entry::Vacant(entry) => {
                            let mut shared = BTreeSet::new();
                            self.collect_owners(definition.rhs.as_view(), &mut shared)?;
                            entry.insert(shared)
                        }
                    };
                    owners.extend(shared.iter().cloned());
                }
                Ok(owners)
            })
            .collect()
    }

    /// Add the owners of the branch tags in `atom` to `owners`.
    fn collect_owners(&self, atom: AtomView<'_>, owners: &mut BTreeSet<SuBitGraph>) -> Result<()> {
        let mut error = None;
        atom.visitor(&mut |view| {
            if error.is_some() {
                return false;
            }
            match self.tag_owner(view) {
                Ok(owner) => {
                    owners.extend(owner);
                    true
                }
                Err(tag_error) => {
                    error = Some(tag_error);
                    false
                }
            }
        });
        error.map_or(Ok(()), Err)
    }

    /// Whether this operation's U acts on a residue row that depends on the
    /// completed zero-mode branches `owners` as the UV counterterm of their
    /// static theory: whether the soft set of the row's outermost branch, the
    /// lines of `loop_lines` outside its support, lies inside this
    /// operation's support. A branch strictly inside another one was taken
    /// over by it, whose zero-mode operation expanded the inner branch's
    /// further soft lines as hard ones. The U then makes every static soft
    /// line of the row massive, so its image of the row is infrared finite on
    /// its own and the masses it introduces there only regulate the static
    /// theory. Otherwise the image keeps massless static lines and the U
    /// introduces the vacuum masses. The own U of a connected host, which
    /// cancels the soft region of its image of the unfrozen integrand and
    /// needs the same masses there, never gets here
    /// ([`Local3DDeformation::OrdinaryOnSoftBranch`]); a U holding part of the
    /// soft set is not constrained in the infrared, since an enclosing host
    /// subtracts the soft region of its image with the image's own masses.
    ///
    /// # Errors
    ///
    /// A row whose outermost branches lie both inside and outside the
    /// operation's support in this sense. Zero-mode hosts form a chain under
    /// inclusion, which leaves a single outermost branch per row.
    pub(super) fn expands_static_cograph(
        &self,
        owners: &BTreeSet<SuBitGraph>,
        loop_lines: &SuBitGraph,
    ) -> Result<bool> {
        let (inside, outside): (Vec<_>, Vec<_>) = owners
            .iter()
            .filter(|owner| {
                !owners
                    .iter()
                    .any(|outer| outer != *owner && outer.includes(*owner))
            })
            .partition(|owner| self.support.includes(&loop_lines.subtract(*owner)));
        if !inside.is_empty() && !outside.is_empty() {
            return Err(eyre!(
                "unsupported zero-mode forest: the U on {} contains the soft sets of the outermost \
                 completed zero-mode branches {} but not those of {} of the same residue row, so \
                 it is neither the UV counterterm of their static theory nor an ordinary U of the \
                 row",
                self.support.string_label(),
                inside.iter().map(|owner| owner.string_label()).join(", "),
                outside.iter().map(|owner| owner.string_label()).join(", "),
            ));
        }
        Ok(!inside.is_empty())
    }

    /// Whether this operation multiplies the tag `view` by its Taylor
    /// parameter, as it does its loop momenta:
    /// - the mass scope or zero-mode temperature of a completed branch iff
    ///   the branch is frozen, being a vertex whose scales follow those
    ///   momenta; the U of a zero-mode branch's own support never sees its
    ///   zero-mode temperature, which only an enclosing zero-mode branch
    ///   takes over;
    /// - the order tag of a frozen zero-mode branch only in an ordinary
    ///   chart, where it compensates the degree of the vertex so that the
    ///   static cograph is expanded at its own degree; a zero-mode chart
    ///   already co-scales the whole vertex;
    /// - the zero-mode temperature of a class iff a zero-mode chart contains
    ///   the class, whose static loop is then part of the hard part it
    ///   expands.
    ///
    /// # Errors
    ///
    /// A tag without a valid owner, a branch partially overlapping the
    /// operation, and a zero-mode temperature in the role of the operation's
    /// own branch.
    pub(super) fn co_scales(&self, view: AtomView<'_>) -> Result<bool> {
        let AtomView::Fun(call) = view else {
            return Ok(false);
        };
        let symbol = call.get_symbol();
        if ![
            *LOCAL_3D_SOFT_SCOPE,
            *ZERO_MODE_TEMPERATURE,
            *ZERO_MODE_ORDER,
        ]
        .contains(&symbol)
        {
            return Ok(false);
        }
        let zero_mode_chart = matches!(self.deformation, Local3DDeformation::ZeroMode { .. });
        match self.tag_role(view)? {
            Some(SoftBranchRole::Frozen) => Ok(symbol != *ZERO_MODE_ORDER || !zero_mode_chart),
            Some(SoftBranchRole::Own) if symbol == *ZERO_MODE_TEMPERATURE => Err(eyre!(
                "internal error: the zero-mode temperature {view} reached the U of its own branch"
            )),
            Some(SoftBranchRole::Own | SoftBranchRole::Disjoint) => Ok(false),
            None => Ok(zero_mode_chart && contains_class(self.graph, self.support, call.get(0))),
        }
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
        let soft_branches = SoftBranches::new(ctx.graph, support, deformation);
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
    /// that were already frozen inside it keep their own owners. A zero-mode
    /// branch also takes over the zero-mode temperatures of the classes and
    /// of the zero-mode branches inside its support, and absorbs the order
    /// tags of the latter: its own order tag now counts the order of the
    /// whole vertex.
    pub(super) fn freeze_soft_branch(
        branch: &DirectResidueBranches,
        graph: &Graph,
        support: &SuBitGraph,
        lmb: &LoopMomentumBasis,
        mass_scope: &Atom,
    ) -> Result<DirectResidueBranches> {
        branch.map_expressions(|atom| {
            Ok(Self::freeze_soft_atom(
                atom, graph, support, lmb, mass_scope,
            ))
        })
    }

    /// [`Self::freeze_soft_branch`] for one atom.
    pub(super) fn freeze_soft_atom(
        atom: &Atom,
        graph: &Graph,
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
        frozen.replace_map(|view, _, output| {
            let AtomView::Fun(call) = view else {
                return;
            };
            let symbol = call.get_symbol();
            if call.get_nargs() != 1
                || ![*ZERO_MODE_TEMPERATURE, *ZERO_MODE_ORDER].contains(&symbol)
            {
                return;
            }
            let tag_owner = call.get(0);
            let inside = match component_owner(tag_owner, graph.n_hedges()) {
                Some(inner) => inner != *support && support.includes(&inner),
                None => {
                    symbol == *ZERO_MODE_TEMPERATURE && contains_class(graph, support, tag_owner)
                }
            };
            if inside {
                **output = if symbol == *ZERO_MODE_TEMPERATURE {
                    zero_mode_temperature(owner.clone())
                } else {
                    Atom::num(1)
                };
            }
        })
    }

    /// Express the soft chart at fixed physical temperature and chemical potential.
    /// For ordinary energy derivatives, N_d(h E; h T, h mu) = h^-d N_d(E; T, mu).
    /// The routed energy is h E(k + q/h, m), so dividing its argument by h
    /// exposes the external-momentum Taylor series. The factor h^-d retains
    /// the derivative's engineering dimension in the existing Laurent counting.
    ///
    /// A soft or zero-mode chart compensates every occupation of its active
    /// support, a union host included. Any deformation also compensates the
    /// occupations retained by a frozen soft branch, keeping their owners,
    /// since that branch is a vertex whose temperature and chemical potential
    /// scale with its loop momenta.
    pub(super) fn compensate_thermal_distributions(
        ctx: &UVCtx<'_>,
        active_subgraph: &SuBitGraph,
        integrand: &Atom,
        deformation: Local3DDeformation,
        soft_branches: SoftBranches<'_>,
    ) -> Result<Atom> {
        // N is the same function for equal statistics and chemical potential.
        // Choose its representative within this active support so coincident
        // energies cancel symbolically, also between the components of a union
        // host: the Landau pairs of the CFF variants then cancel in each
        // Taylor coefficient. Ownership of the branch lives in the energy's
        // owner, not in this representative; disjoint components and their
        // occupations remain unchanged.
        let mut representatives = BTreeMap::new();
        let owners = if deformation.is_soft() {
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
