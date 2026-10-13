use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Neg,
    sync::Arc,
    time::Instant,
};

use color_eyre::Result;
use eyre::eyre;
use idenso::color::{ColorSimplifier, ColorSimplifySettings};
use linnet::half_edge::{
    involution::{EdgeIndex, Flow},
    subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate},
    domains::rational::Rational,
    function,
    id::Replacement,
    symbol,
};
use three_dimensional_reps::CffGenerationContext;

use crate::{
    cff::{
        CutCFFIndex,
        expression::OrientationID,
        surface::{GammaLoopLinearEnergyExpr, LinearEnergyExpr},
    },
    graph::{ExactUvSubLmbFrame, Graph, LMBext, LoopMomentumBasis, cuts::CutSet},
    integrands::process::param_builder::FnMapEntry,
    utils::{GS, W_},
    uv::{
        Integrands, UVgenerationSettings, UltravioletGraph,
        approx::{
            ForestNodeLike, UVCtx,
            integrated::IntegratedCts,
            local_3d::{
                FrozenActiveCt, Local3DLoopRescaling, Local3DProjectionPath, Localizer,
                paths_with_step,
            },
            local_4d::FourDTerm,
            projected_4d::Local4dProjectionContext,
        },
        hedge_poset::OwnedForestNode,
        marker::{UvMarker, UvOperation},
        spinney::Spinney,
        uv_graph::UVE,
    },
};

use super::{
    branches::{DirectResidueBranches, DirectResidueKey},
    kernel::{DirectCoordinateFrame, LOCAL_3D_MASS_SCOPE, apply_taylor, coordinate_lmb, start},
};

fn extend_coordinate_frames(
    prior_frames: Vec<DirectCoordinateFrame>,
    frame: DirectCoordinateFrame,
) -> Result<Vec<DirectCoordinateFrame>> {
    let mut next = Vec::with_capacity(prior_frames.len() + 1);
    for prior in prior_frames {
        if frame.active_subgraph.intersects(&prior.active_subgraph) {
            if !frame.active_subgraph.includes(&prior.active_subgraph) {
                return Err(eyre!(
                    "a direct local-3D coordinate frame partially overlaps a prior forest frame"
                ));
            }
        } else {
            next.push(prior);
        }
    }
    next.push(frame);
    Ok(next)
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct DirectSector {
    pub(crate) active_subgraph: SuBitGraph,
    /// Loop coordinates in which each connected part of `active` was Taylor
    /// expanded. Descendants retain disjoint frames separately and combine
    /// them only when one enclosing Taylor operation contains both parts.
    pub(crate) coordinate_frames: Vec<DirectCoordinateFrame>,
    pub(crate) projection_paths: Vec<Local3DProjectionPath>,
    pub(crate) active: DirectResidueBranches,
    pub(crate) frozen_integrands: Integrands,
}

impl DirectSector {
    pub(crate) fn combine(&self) -> Result<DirectResidueBranches> {
        // Scope is needed by later forest operations only. Materialize the physical
        // vacuum mass on this output copy while retaining owners in the stored sector.
        let physical = self.active.map_expressions(|atom| {
            Ok(atom
                .replace(function!(GS.m_uv_vacuum, W_.x_))
                .with(GS.m_uv_vacuum))
        })?;
        physical.zip_mul_unmapped(&self.frozen_integrands)
    }
}

impl Neg for DirectSector {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            active_subgraph: self.active_subgraph,
            coordinate_frames: self.coordinate_frames,
            projection_paths: self.projection_paths,
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
        // Bridge this source's CFF energy-factor convention once at the root.
        // Forest Taylor subtraction signs are applied separately below.
        let production_prefactor = Atom::num(cff.production_prefactor_factor());
        let exact_orientations = localizer.orientation.exact_orientations()?;
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
                    GS.retain_esurfaces(&stored.expression) * &production_prefactor;
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

    pub(crate) fn projection_paths(&self) -> Vec<Local3DProjectionPath> {
        match self {
            Self::Root(_) => Vec::new(),
            Self::Sectors(sectors) => {
                sectors
                    .iter()
                    // Retain zero sectors for active-mask replay, but do not publish
                    // a projection history for a branch that contributes no integrand.
                    .filter(|sector| {
                        sector.active.iter_keys().any(|(_, integrands)| {
                            integrands.iter().any(|(residue, atom)| {
                                !atom.is_zero()
                                    && sector.frozen_integrands.iter().any(
                                        |(frozen_residue, frozen)| {
                                            residue == frozen_residue && !frozen.is_zero()
                                        },
                                    )
                            })
                        })
                    })
                    .flat_map(|sector| sector.projection_paths.iter().cloned())
                    .collect()
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
                            coordinate_frames: sector.coordinate_frames.clone(),
                            projection_paths: sector.projection_paths.clone(),
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
    /// The contour on which this replay operates. Production forests use the
    /// whole graph; a completed proper O source replays within its own frame
    /// before its coefficient is attached to the live surrounding contour.
    source_frame: DirectCoordinateFrame,
}

/// One component's ordered projector branches and physically linked fixed
/// point, before attaching any surrounding source or completed coefficient.
struct DirectOnShellPreparation {
    branches: [Atom; 2],
    owner: EdgeIndex,
    scope_edge: EdgeIndex,
    mass: Atom,
    mass_scope: Atom,
    physical_masses: BTreeSet<Atom>,
    fixed_coordinates: Vec<(EdgeIndex, i32)>,
}

impl DirectOnShellPreparation {
    fn new<S: ForestNodeLike>(
        graph: &Graph,
        current: &S,
        given: &S,
        lmb: &LoopMomentumBasis,
    ) -> Result<Self> {
        let numerator = graph
            .numerator(&current.reduced_subgraph(given), given.subgraph())
            .get_single_atom()
            .expect("the component numerator is available")
            .simplify_color_with(ColorSimplifySettings {
                simplify_non_color: false,
                ..Default::default()
            });
        let branches = graph.on_shell_branches(&numerator, current.subgraph(), lmb)?;
        let crown = graph.dummy_stripped_external_flows_of(current.subgraph());
        let incoming = crown
            .included_iter()
            .find(|hedge| graph.flow(*hedge) == Flow::Sink)
            .ok_or_else(|| eyre!("direct OS source requires an incoming quark boundary"))?;
        let owner = graph[&incoming];
        let scope_edge = *lmb
            .loop_edges
            .first()
            .ok_or_else(|| eyre!("OS source has no loop coordinate"))?;
        let owner_signature = &graph.loop_momentum_basis.edge_signatures[owner];
        let fixed_coordinates = lmb
            .ext_edges
            .iter()
            .filter_map(|edge| {
                let signature = &graph.loop_momentum_basis.edge_signatures[*edge];
                signature
                    .equality_up_to_sign(owner_signature)
                    .then_some((*edge, if signature == owner_signature { 1 } else { -1 }))
            })
            .collect::<Vec<_>>();
        if fixed_coordinates.is_empty() {
            return Err(eyre!(
                "OS source has no exact external coordinate for its incoming point"
            ));
        }
        Ok(Self {
            branches,
            owner,
            scope_edge,
            mass: graph[owner].mass_atom(),
            mass_scope: function!(*LOCAL_3D_MASS_SCOPE, usize::from(scope_edge)),
            physical_masses: graph
                .iter_edges_of(current.subgraph())
                .filter(|(pair, _, _)| pair.is_paired())
                .map(|(_, edge, _)| graph[edge].mass_atom())
                .filter(|mass| !mass.is_zero())
                .collect(),
            fixed_coordinates,
        })
    }

    fn scope_masses(&self, atom: &Atom) -> Atom {
        atom.replace_multiple(self.physical_masses.iter().map(|mass| {
            Replacement::new(mass.to_pattern(), (mass * &self.mass_scope).to_pattern())
        }))
    }

    fn point(&self, sigma: i32) -> Atom {
        function!(
            GS.os_mom,
            usize::from(self.owner),
            Atom::num(sigma) * &self.mass * &self.mass_scope
        )
    }

    /// A complete direct residue map samples original numerator occurrences.
    /// Equal routed momenta on a raised channel can have different contact
    /// samples, so restore the NEW numerator's owner before that map acts.
    /// The complete-domain source binds its external point after mapping.
    fn original_source_numerator(&self, numerator: &Atom) -> Atom {
        let original = numerator.replace_map(|view, _, output| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() == 0 {
                return;
            }
            let Some((
                owner,
                crate::utils::symbols::UvMomentumProvenanceRole::PhysicalSourceFixed,
                _,
            )) = GS.uv_momentum_provenance_data(momentum.get(0))
            else {
                return;
            };
            let indices = momentum
                .iter()
                .skip(1)
                .map(|index| index.to_owned())
                .collect::<Vec<_>>();
            **output = GS.indexed_momentum(&function!(GS.emr_mom, usize::from(owner)), &indices);
        });
        self.scope_masses(&original)
    }
}

impl<'a> Direct3dApproximation<'a> {
    pub(crate) fn new(
        localizer: Localizer<'a>,
        graph: &'a mut Graph,
        settings: &'a UVgenerationSettings,
    ) -> Self {
        let source_frame = DirectCoordinateFrame {
            active_subgraph: graph.full_filter(),
            lmb: graph.loop_momentum_basis.clone(),
        };
        Self {
            localizer,
            graph,
            settings,
            source_frame,
        }
    }

    /// O fixes only the current component's external data. A proper primitive
    /// component owns its original contour and a separate live cograph. When
    /// the current component spans the graph, its complete incoming residue
    /// family is the source, including every completed child coefficient.
    pub(super) fn on_shell_source<S: ForestNodeLike>(
        &mut self,
        current: &S,
        given: &S,
        lmb: &LoopMomentumBasis,
        active: &DirectResidueBranches,
        active_subgraph: &SuBitGraph,
    ) -> Result<DirectResidueBranches> {
        let complete_domain = self.graph.paired_edges(current.subgraph())
            == self.graph.paired_edges(&self.source_frame.active_subgraph);
        if current.dod() != 1 || (!complete_domain && !given.subgraph().is_empty()) {
            return Err(eyre!(
                "direct OS requires a complete current domain or an original component source"
            ));
        }
        if !complete_domain {
            if self.graph.n_loops(current.subgraph()) != 1 {
                return self.proper_on_shell_source(current, given, lmb, None);
            }
            let coefficient = self.on_shell_component(current, given, lmb)?;
            return self.attach_on_shell_components(current.subgraph(), &coefficient, &Atom::one());
        }
        let prepared = DirectOnShellPreparation::new(self.graph, current, given, lmb)?;
        // Scalar definitions may themselves contain parent external data.
        // Resolve only this scalar layer, preserving shared tensor bodies
        // and the frozen localization kernels held outside `active`.
        let active = active.resolved_scalars()?;
        let ctx = UVCtx::new(self.graph, self.settings);
        let mut completed = Vec::with_capacity(2);
        for (sigma, numerator) in [1, -1].into_iter().zip(&prepared.branches) {
            let numerator = prepared.original_source_numerator(numerator);
            let branch = active.multiply_key_mapped(
                self.localizer.orientation,
                self.graph,
                &numerator,
                DirectResidueBranches::numerator_scope(),
            )?;
            let point = prepared.point(sigma);
            completed.push(branch.map_expressions(|atom| {
                let started = start(
                    &ctx,
                    current,
                    atom,
                    &Atom::one(),
                    Some(active_subgraph),
                    lmb,
                )?;
                let routed = started.replace_multiple(self.graph.uv_spatial_wrapped_replacement(
                    current.subgraph(),
                    lmb,
                    &[W_.x___],
                ));
                let evaluated = routed.replace_map(|view, _, output| {
                    let AtomView::Fun(momentum) = view else {
                        return;
                    };
                    if ![GS.emr_mom, GS.emr_vec].contains(&momentum.get_symbol())
                        || momentum.get_nargs() == 0
                    {
                        return;
                    }
                    let Ok(edge) = usize::try_from(momentum.get(0)) else {
                        return;
                    };
                    let Some((_, sign)) = prepared
                        .fixed_coordinates
                        .iter()
                        .find(|(coordinate, _)| usize::from(*coordinate) == edge)
                    else {
                        return;
                    };
                    **output = if momentum.get_symbol() == GS.emr_vec {
                        Atom::Zero
                    } else {
                        let indices = momentum
                            .iter()
                            .skip(1)
                            .map(|index| index.to_owned())
                            .collect::<Vec<_>>();
                        Atom::num(*sign) * GS.indexed_momentum(&point, &indices)
                    };
                });
                Ok(GS.resolve_os_momenta(&evaluated))
            })?);
        }
        completed[0].zip_add(&completed[1])
    }

    /// Integrate only the original physical O branches of one component. The
    /// ordinary U operation continues to act on its hosted direct CFF family.
    fn on_shell_component<S: ForestNodeLike>(
        &mut self,
        current: &S,
        given: &S,
        lmb: &LoopMomentumBasis,
    ) -> Result<Integrands> {
        let prepared = DirectOnShellPreparation::new(self.graph, current, given, lmb)?;
        let crown = self
            .graph
            .dummy_stripped_external_flows_of(current.subgraph());
        // This parses only the original propagator product. No U, S, or
        // recursive Local4dCts object is used to create this direct contour.
        let mut source =
            FourDTerm::from_view(self.graph.denominator(current.subgraph(), |_| -1).as_view())?;
        let [source] = source.as_mut_slice() else {
            return Err(eyre!(
                "an original component must have one denominator product"
            ));
        };
        for denominator in &mut source.denominators {
            denominator.mass_squared = prepared.scope_masses(&denominator.mass_squared);
            denominator.full_expr = prepared.scope_masses(&denominator.full_expr);
        }
        let compact_energies = |atom: &Atom| {
            atom.replace_map(|view, _, output| {
                let AtomView::Pow(power) = view else {
                    return;
                };
                let Ok(exponent) = Rational::try_from(power.get_exp()) else {
                    return;
                };
                if exponent.denominator() != 2 {
                    return;
                }
                // These are the positive source pole energies. Keeping their
                // explicit invariant in the established compact dialect lets a
                // subsequent direct U apply its mass rearrangement to this loop.
                if !power.get_base().contains_symbol(GS.emr_mom)
                    && !power.get_base().contains_symbol(GS.emr_vec)
                    && !power.get_base().contains_symbol(*LOCAL_3D_MASS_SCOPE)
                {
                    return;
                }
                **output = function!(
                    GS.on_shell_energy,
                    usize::from(prepared.scope_edge),
                    power.get_base()
                )
                .pow(Atom::num(exponent.numerator()));
            })
        };
        let mut options = self.localizer.orientation.cff_options()?;
        options.cff_generation_context = CffGenerationContext::EmbeddedCffFactor;
        let child_cutset = CutSet::empty(self.graph.n_hedges());
        let mut context = Local4dProjectionContext::default();
        let mut coefficients = Vec::new();
        for (sigma, numerator) in [1, -1].into_iter().zip(&prepared.branches) {
            let numerator = prepared.scope_masses(&GS.erase_uv_momentum_provenance(numerator));
            let point = prepared.point(sigma);
            let fixed = prepared
                .fixed_coordinates
                .iter()
                .map(|(edge, sign)| (*edge, Atom::num(*sign) * &point))
                .collect();
            let (cff, _) = self.graph.cff_from_4d_denominators_in_uv_sub_lmb(
                &source.denominators,
                self.graph.paired_edges(current.subgraph()),
                crown.included_iter(),
                lmb,
                ExactUvSubLmbFrame::RetainedPhysicalCrown,
                &child_cutset,
                &options,
                &numerator,
                Some(&mut context),
                &[],
                &fixed,
            )?;
            let prefactor = Atom::num(cff.production_prefactor_factor());
            let scope = DirectResidueBranches::numerator_scope().1;
            let mut coefficient = Atom::Zero;
            let mut definitions = Vec::new();
            for (cut, term) in cff.terms {
                if cut != CutCFFIndex::new_all_none() {
                    return Err(eyre!("an OS child contour must be uncut"));
                }
                let (body, parameters, rows) = term.prepare_exact_source_numerator(&scope)?;
                let lhs = symbol!("gammalooprs::uv::numerator_family")
                    .call_args(std::iter::once(scope.clone()).chain(parameters.iter().cloned()));
                let entry = Arc::new(FnMapEntry {
                    lhs,
                    rhs: compact_energies(&body),
                    args: parameters
                        .iter()
                        .cloned()
                        .map(Indeterminate::try_from)
                        .collect::<std::result::Result<Vec<_>, _>>()
                        .map_err(|error| eyre!(error))?,
                    tags: vec![scope.clone()],
                    inlining: symbolica::evaluate::InliningPolicy::Always,
                    is_alias: false,
                });
                for (orientation, row) in term.orientations.iter().zip(rows) {
                    let call = entry.lhs.replace_multiple(parameters.iter().zip(row).map(
                        |(parameter, value)| {
                            Replacement::new(
                                parameter.to_pattern(),
                                compact_energies(&value).to_pattern(),
                            )
                        },
                    ));
                    coefficient += compact_energies(&orientation.expression) * &prefactor * call;
                }
                definitions.push(entry);
            }
            coefficients.push(
                Integrands::from_iter([(CutCFFIndex::new_all_none(), coefficient)])
                    .with_numerators(definitions)?,
            );
        }
        let mut coefficients = coefficients.into_iter();
        coefficients
            .next()
            .expect("two OS branches")
            .zip_add(coefficients)
    }

    /// Attach independent O coefficients to their live quotient exactly once.
    fn attach_on_shell_components(
        &mut self,
        contracted: &SuBitGraph,
        coefficient: &Integrands,
        finite_coefficient: &Atom,
    ) -> Result<DirectResidueBranches> {
        let remaining = self.source_frame.active_subgraph.subtract(contracted);
        let cograph_numerator = self
            .graph
            .numerator(&remaining, contracted)
            .get_single_atom()
            .expect("the cograph numerator is available")
            * self.graph.global_atom();
        let analysis = coefficient
            .numerators()
            .iter()
            .map(|entry| &entry.rhs * &cograph_numerator * finite_coefficient)
            .collect::<Vec<_>>();
        let outer = self
            .source_in_frame(contracted, &analysis)?
            .multiply_key_mapped(
                self.localizer.orientation,
                self.graph,
                finite_coefficient,
                DirectResidueBranches::numerator_scope(),
            )?;
        let scalar = coefficient.iter().next().expect("uncut OS coefficient").1;
        let broadcast = Integrands::from_iter(
            outer
                .identity_integrands()
                .iter()
                .map(|(key, _)| (*key, scalar.clone())),
        )
        .with_numerators(coefficient.numerators().iter().cloned())?;
        let tree = GS.wrap_tree_denoms(
            self.graph.denominator(
                &self
                    .graph
                    .tree_edges
                    .intersection(&self.source_frame.active_subgraph)
                    .subtract(&self.graph.initial_state_cut),
                |_| -1,
            ),
        );
        Ok(outer.zip_mul_unmapped(&broadcast)?.map(|atom| atom * &tree))
    }

    /// Construct a denominator-only source in the same explicit domain as
    /// the replay. Its full physical numerator bounds and certifies sampling,
    /// but each forest operation remains the sole owner of its numerator.
    fn source_in_frame(
        &mut self,
        contracted: &SuBitGraph,
        analysis: &[Atom],
    ) -> Result<DirectResidueBranches> {
        if self.graph.paired_edges(&self.source_frame.active_subgraph)
            == self.graph.paired_edges(&self.graph.full_filter())
        {
            let source = self.localizer.with_independent_source_sum().projected_cff(
                self.graph,
                contracted,
                analysis,
                CffGenerationContext::EmbeddedCffFactor,
            )?;
            return DirectResidueBranches::from_transient(&source);
        }
        if !self.source_frame.active_subgraph.includes(contracted) || analysis.is_empty() {
            return Err(eyre!(
                "a component source requires contained owners and a numerator witness"
            ));
        }
        let remaining = self.source_frame.active_subgraph.subtract(contracted);
        let crown = self
            .graph
            .dummy_stripped_external_flows_of(&self.source_frame.active_subgraph);
        let basis = if contracted.is_empty() {
            self.source_frame.lmb.clone()
        } else {
            let contracted = InternalSubGraph::try_new(contracted.clone(), self.graph.as_ref())
                .ok_or_else(|| eyre!("a component source can contract only internal prefixes"))?;
            self.graph.shrunken_sub_lmb(
                &self.source_frame.active_subgraph,
                &contracted,
                crown.clone(),
                Some(&self.source_frame.lmb),
            )?
        };
        let source = FourDTerm::from_view(self.graph.denominator(&remaining, |_| -1).as_view())?;
        let [source] = source.as_slice() else {
            return Err(eyre!(
                "an original component source requires one propagator product"
            ));
        };
        // Distinguish alternative coefficient bodies in the capacity witness.
        // This prevents cancellations between them from lowering the contour
        // rank. These inert weights never enter the emitted denominator CFF.
        let scope = DirectResidueBranches::numerator_scope().1;
        let numerator = analysis
            .iter()
            .enumerate()
            .fold(Atom::Zero, |sum, (index, atom)| {
                sum + function!(
                    symbol!("gammalooprs::uv::source_capacity"; Scalar),
                    &scope,
                    index
                ) * atom
            });
        let replacements = self.graph.uv_spatial_wrapped_replacement(
            &self.source_frame.active_subgraph,
            &basis,
            &[W_.x___],
        );
        let mut poles = BTreeMap::new();
        for denominator in &source.denominators {
            let invariant = (1..=3)
                .fold(denominator.mass_squared.clone(), |sum, component| {
                    sum + GS
                        .indexed_momentum(&denominator.momentum, &[GS.cind(component)])
                        .pow(2)
                })
                .replace_multiple(&replacements)
                .expand();
            poles.entry(invariant).or_insert(denominator.source_edge);
        }
        fn quadratic_degree(view: AtomView<'_>) -> Option<usize> {
            let degree = match view {
                AtomView::Num(_) => 0,
                AtomView::Var(_) => 1,
                AtomView::Fun(function) if function.get_symbol() == GS.emr_mom => 1,
                AtomView::Add(sum) => sum
                    .iter()
                    .map(quadratic_degree)
                    .collect::<Option<Vec<_>>>()?
                    .into_iter()
                    .max()
                    .unwrap_or(0),
                AtomView::Mul(product) => product.iter().try_fold(0usize, |degree, factor| {
                    Some(degree + quadratic_degree(factor)?)
                })?,
                AtomView::Pow(power) => quadratic_degree(power.get_base())?
                    .checked_mul(usize::try_from(power.get_exp()).ok()?)?,
                _ => return None,
            };
            (degree <= 2).then_some(degree)
        }
        let compact = |atom: &Atom| {
            atom.replace_map(|view, _, output| {
                let AtomView::Pow(power) = view else {
                    return;
                };
                let Ok(exponent) = Rational::try_from(power.get_exp()) else {
                    return;
                };
                // Only a certified small quadratic can be a source pole
                // invariant. Do not expand a causal surface or a graph sum
                // while searching for an energy alias.
                if quadratic_degree(power.get_base()).is_none() {
                    return;
                }
                if let Some(edge) = poles.get(&power.get_base().expand()) {
                    **output = GS.ose(*edge).pow(Atom::num(exponent * Rational::from(2)));
                }
            })
        };
        let mut options = self.localizer.orientation.cff_options()?;
        options.cff_generation_context = CffGenerationContext::EmbeddedCffFactor;
        let cutset = CutSet::empty(self.graph.n_hedges());
        let (cff, _) = self.graph.cff_from_4d_denominators_in_uv_sub_lmb(
            &source.denominators,
            self.graph.paired_edges(&self.source_frame.active_subgraph),
            crown.included_iter(),
            &basis,
            ExactUvSubLmbFrame::RetainedPhysicalCrown,
            &cutset,
            &options,
            &numerator,
            None,
            &[],
            &BTreeMap::new(),
        )?;
        let prefactor = Atom::num(cff.production_prefactor_factor());
        let host = self
            .localizer
            .orientation
            .orientation_ids()?
            .into_iter()
            .next()
            .ok_or_else(|| eyre!("a component source requires a permitted selector host"))?;
        let mut branches = Vec::new();
        for (cut, term) in cff.terms {
            if cut != CutCFFIndex::new_all_none() {
                return Err(eyre!("the internal component replay must be uncut"));
            }
            let rows = term.certified_original_source_energy_rows(
                &numerator,
                self.graph.underlying.n_edges(),
            )?;
            for (orientation, row) in term.orientations.into_iter().zip(rows) {
                let map = row
                    .iter()
                    .map(|energy| LinearEnergyExpr::from_atom_gs(&compact(energy)))
                    .collect::<Result<Vec<_>>>()?;
                branches.push((
                    DirectResidueKey::source(host, map),
                    Integrands::from_iter([(
                        cut,
                        GS.retain_esurfaces(&compact(&orientation.expression)) * &prefactor,
                    )]),
                ));
            }
        }
        DirectResidueBranches::from_keyed(branches)
    }

    /// Commute O through one disjoint completed primitive operation. Its U
    /// remains on the original hosted family; only the independent O factors
    /// own separate contours. No already completed coefficient is discarded.
    #[allow(clippy::too_many_arguments)]
    fn disjoint_on_shell_source<S: ForestNodeLike>(
        &mut self,
        current: &S,
        given: &S,
        lmb: &LoopMomentumBasis,
        prior_active: &SuBitGraph,
        frames: &[DirectCoordinateFrame],
        paths: &[Local3DProjectionPath],
        frozen: &Integrands,
    ) -> Result<DirectResidueBranches> {
        if prior_active.intersects(current.subgraph())
            || !given.subgraph().is_empty()
            || self.graph.n_loops(current.subgraph()) != 1
            || current.dod() != 1
        {
            return Err(eyre!(
                "disjoint direct OS replay requires independent primitive components"
            ));
        }
        let prior = self.primitive_local_prefix(prior_active, frames, paths, frozen)?;
        let root = OwnedForestNode {
            spinney: Spinney::empty(self.graph),
            topo_order: 0,
        };
        let right = self.on_shell_component(current, given, lmb)?;
        let right_source =
            self.attach_on_shell_components(current.subgraph(), &right, &Atom::one())?;
        let both = if prior.renormalization_scheme() == crate::uv::ApproximationType::OS {
            let left = self.on_shell_component(&prior, &root, prior.lmb())?;
            Some(self.attach_on_shell_components(
                &prior_active.union(current.subgraph()),
                &left.zip_mul(&right)?,
                &Atom::one(),
            )?)
        } else {
            None
        };
        let replay = -apply_taylor(
            &UVCtx::new(self.graph, self.settings),
            self.localizer.orientation,
            &prior,
            &root,
            Some(prior_active.clone()),
            prior.lmb(),
            &right_source,
            both.as_ref(),
        )?;
        Ok(replay.map(|atom| {
            UvMarker::new(self.settings).apply(
                UvOperation::Approx,
                prior.subgraph(),
                root.subgraph(),
                atom,
            )
        }))
    }

    fn primitive_local_prefix(
        &self,
        prior_active: &SuBitGraph,
        frames: &[DirectCoordinateFrame],
        paths: &[Local3DProjectionPath],
        frozen: &Integrands,
    ) -> Result<OwnedForestNode> {
        let ([frame], [path]) = (frames, paths) else {
            return Err(eyre!(
                "direct component replay requires one certified primitive prefix"
            ));
        };
        let [step] = path.steps.as_slice() else {
            return Err(eyre!(
                "direct component replay requires one primitive operation"
            ));
        };
        let empty: SuBitGraph = self.graph.empty_subgraph();
        if &frame.active_subgraph != prior_active
            || step.current_component != prior_active.string_label()
            || step.active_subgraph != prior_active.string_label()
            || step.rescaled_subgraph != prior_active.string_label()
            || step.given_component != empty.string_label()
            || step.rescaling != Local3DLoopRescaling::FullSubgraph
            || self.graph.n_loops(prior_active) != 1
            || !frozen.iter().all(|(_, atom)| atom.is_one())
        {
            return Err(eyre!(
                "direct component replay requires a primitive local prefix; finite seeds retain their explicit integrated owner"
            ));
        }
        if step.route_loop_edges
            != frame
                .lmb
                .loop_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .collect::<Vec<_>>()
            || step.route_external_edges
                != frame
                    .lmb
                    .ext_edges
                    .iter()
                    .map(|edge| usize::from(*edge))
                    .collect::<Vec<_>>()
            || step.route_signatures.len() != self.graph.paired_edges(prior_active).len()
            || step.route_signatures.iter().any(|record| {
                let signature = &frame.lmb.edge_signatures[EdgeIndex(record.edge)];
                record.loop_signature != signature.internal.to_string()
                    || record.external_signature != signature.external.to_string()
            })
        {
            return Err(eyre!(
                "the direct prefix changed its recorded momentum frame"
            ));
        }
        let mut spinney = Spinney::with_scheme(
            InternalSubGraph::try_new(prior_active.clone(), self.graph.as_ref())
                .ok_or_else(|| eyre!("the direct prefix is not an internal subgraph"))?,
            self.graph,
            &frame.lmb,
            step.scheme,
            step.dod,
        )?
        .ok_or_else(|| eyre!("the recorded direct prefix is not divergent"))?;
        spinney.lmb = frame.lmb.clone();
        Ok(OwnedForestNode {
            spinney,
            topo_order: step.topo_order,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn nested_on_shell_source<S: ForestNodeLike>(
        &mut self,
        current: &S,
        given: &S,
        lmb: &LoopMomentumBasis,
        prior_active: &SuBitGraph,
        frames: &[DirectCoordinateFrame],
        paths: &[Local3DProjectionPath],
        frozen: &Integrands,
    ) -> Result<DirectResidueBranches> {
        if !current.subgraph().includes(prior_active)
            || given.subgraph() != prior_active
            || current.dod() != 1
        {
            return Err(eyre!(
                "proper OS replay requires a contained, complete recorded prefix"
            ));
        }
        let prior = self.primitive_local_prefix(prior_active, frames, paths, frozen)?;
        self.proper_on_shell_source(current, given, lmb, Some(&prior))
    }

    /// Reuse the same original component contour for a bare component and
    /// for a complete recorded local prefix. The live outer graph is attached
    /// only after the component's fixed-point operation has finished.
    fn proper_on_shell_source<S: ForestNodeLike>(
        &mut self,
        current: &S,
        given: &S,
        lmb: &LoopMomentumBasis,
        prior: Option<&OwnedForestNode>,
    ) -> Result<DirectResidueBranches> {
        let root = OwnedForestNode {
            spinney: Spinney::empty(self.graph),
            topo_order: 0,
        };
        let numerator = self
            .graph
            .numerator(current.subgraph(), root.subgraph())
            .get_single_atom()
            .expect("the original parent numerator is available");
        let coefficient = {
            let mut component = Direct3dApproximation {
                localizer: self.localizer,
                graph: &mut *self.graph,
                settings: self.settings,
                source_frame: DirectCoordinateFrame {
                    active_subgraph: current.subgraph().clone(),
                    lmb: lmb.clone(),
                },
            };
            let source =
                Direct3dCts::Root(component.source_in_frame(root.subgraph(), &[numerator])?);
            let active = if let Some(prior) = prior {
                let replay = component.run_local(&source, prior, &root, prior, &root)?;
                let [replay] = replay.sectors()? else {
                    return Err(eyre!(
                        "one primitive component replay must retain one local sector"
                    ));
                };
                replay.active.clone()
            } else {
                source.branches()?
            };
            component
                .on_shell_source(current, given, lmb, &active, current.subgraph())?
                .materialize(false)?
        };
        self.attach_on_shell_components(current.subgraph(), &coefficient, &Atom::one())
    }

    /// Convert the projection helper's transient representation into the
    /// direct lane's explicit full-key type. A source map remains attached to
    /// its branch and therefore governs every later factor in that branch.
    fn localize_integrated<S: ForestNodeLike>(
        &mut self,
        expr: &Atom,
        integrated_node: &S,
    ) -> Result<(DirectResidueBranches, Integrands)> {
        if self.graph.paired_edges(&self.source_frame.active_subgraph)
            != self.graph.paired_edges(&self.graph.full_filter())
        {
            if !self
                .source_frame
                .active_subgraph
                .includes(integrated_node.subgraph())
            {
                return Err(eyre!(
                    "a component finite source cannot contract a foreign owner"
                ));
            }
            let remaining = self
                .source_frame
                .active_subgraph
                .subtract(integrated_node.subgraph());
            let outside = self
                .graph
                .numerator(&remaining, integrated_node.subgraph())
                .get_single_atom()
                .expect("the reduced component numerator is available");
            let source = self.source_in_frame(integrated_node.subgraph(), &[expr * outside])?;
            let active = source.multiply_key_mapped(
                self.localizer.orientation,
                self.graph,
                expr,
                DirectResidueBranches::numerator_scope(),
            )?;
            let tree = GS.wrap_tree_denoms(
                self.graph.denominator(
                    &self
                        .graph
                        .tree_edges
                        .intersection(&self.source_frame.active_subgraph)
                        .subtract(&self.graph.initial_state_cut),
                    |_| -1,
                ),
            );
            let frozen = active
                .identity_integrands()
                .map(|_| GS.localizing_integrand(integrated_node.lmb()));
            return Ok((active.map(|atom| atom * &tree), frozen));
        }
        let FrozenActiveCt {
            active,
            frozen_integrands,
        } = self.localizer.localize(expr, self.graph, integrated_node)?;
        let branches = DirectResidueBranches::from_transient(&active)?;
        Ok((branches, frozen_integrands))
    }

    #[allow(clippy::too_many_arguments)]
    fn nested_finite_on_shell_source<S: ForestNodeLike, I: ForestNodeLike>(
        &mut self,
        current: &S,
        given: &S,
        integrated_node: &I,
        finite_coefficient: &Atom,
        canonical_lmb: &LoopMomentumBasis,
        coordinate_lmb: &LoopMomentumBasis,
        active_subgraph: &SuBitGraph,
    ) -> Result<DirectResidueBranches> {
        let coefficient = {
            let mut component = Direct3dApproximation {
                localizer: self.localizer,
                graph: &mut *self.graph,
                settings: self.settings,
                source_frame: DirectCoordinateFrame {
                    active_subgraph: current.subgraph().clone(),
                    // This frame owns the full component. Its source builder
                    // contracts the given prefix once; O itself uses the
                    // resulting active quotient coordinates below.
                    lmb: canonical_lmb.clone(),
                },
            };
            let (active, _) = component.localize_integrated(finite_coefficient, integrated_node)?;
            component
                .on_shell_source(current, given, coordinate_lmb, &active, active_subgraph)?
                .materialize(false)?
        };
        // The caller retains the original normalized child kernel outside
        // the operator. This coefficient contains only its analytic seed and
        // the unexpanded active parent contour.
        self.attach_on_shell_components(current.subgraph(), &coefficient, &Atom::one())
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
        // Localize the integrated source before reading the local graph, as
        // localization registers its CFF surfaces. The two Taylor branches
        // are independent; each replay method supplies its own minus and marker.
        // A sequential zero addback is omitted, while an explicit integrated
        // replay keeps its typed zero sector for disconnected descendants.
        let finite_counterterm = integrated.physical_finite_counterterm_atom();
        let integrated = (!finite_counterterm.is_zero())
            .then(|| {
                self.run_integrated(
                    &[(&finite_counterterm, given)],
                    given,
                    current,
                    given,
                    marker_current,
                    marker_given,
                )
            })
            .transpose()?;
        let local = self.run_local(local, current, given, marker_current, marker_given)?;
        if let Some(integrated) = integrated {
            let mut sectors = local.sectors()?.to_vec();
            sectors.extend(integrated.sectors()?.iter().cloned());
            Direct3dCts::from_sectors(sectors)
        } else {
            Ok(local)
        }
    }

    pub(crate) fn run_local<S: ForestNodeLike, M: ForestNodeLike>(
        &mut self,
        local: &Direct3dCts,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Direct3dCts> {
        let started = Instant::now();
        let reduced_subgraph = current.reduced_subgraph(given);
        let sectors = match local {
            Direct3dCts::Root(branches) => vec![(
                self.graph.empty_subgraph(),
                Vec::new(),
                Vec::new(),
                branches,
                branches.identity_integrands(),
            )],
            Direct3dCts::Sectors(sectors) => sectors
                .iter()
                .map(|sector| {
                    (
                        sector.active_subgraph.clone(),
                        sector.coordinate_frames.clone(),
                        sector.projection_paths.clone(),
                        &sector.active,
                        sector.frozen_integrands.clone(),
                    )
                })
                .collect(),
        };
        let next = sectors
            .into_iter()
            .map(
                |(prior_active_subgraph, prior_frames, prior_paths, active, frozen_integrands)| {
                    let active_subgraph = prior_active_subgraph.union(&reduced_subgraph);
                    // Retain the complete active mask for descendants, but rescale
                    // only the component-local part covered by this path.
                    let rescaled_subgraph = active_subgraph.intersection(current.subgraph());
                    let (canonical_lmb, coordinate_lmb) = coordinate_lmb(
                        &UVCtx::new(self.graph, self.settings),
                        current,
                        given,
                        Some(&prior_active_subgraph),
                        &prior_frames,
                        &rescaled_subgraph,
                    )?;
                    let on_shell = if current.renormalization_scheme()
                        == crate::uv::ApproximationType::OS
                    {
                        Some(
                            if !prior_active_subgraph.is_empty()
                                && self.graph.paired_edges(current.subgraph())
                                    != self.graph.paired_edges(&self.source_frame.active_subgraph)
                            {
                                if current.subgraph().includes(&prior_active_subgraph) {
                                    self.nested_on_shell_source(
                                        current,
                                        given,
                                        &coordinate_lmb,
                                        &prior_active_subgraph,
                                        &prior_frames,
                                        &prior_paths,
                                        &frozen_integrands,
                                    )?
                                } else {
                                    self.disjoint_on_shell_source(
                                        current,
                                        given,
                                        &coordinate_lmb,
                                        &prior_active_subgraph,
                                        &prior_frames,
                                        &prior_paths,
                                        &frozen_integrands,
                                    )?
                                }
                            } else {
                                self.on_shell_source(
                                    current,
                                    given,
                                    &coordinate_lmb,
                                    active,
                                    &rescaled_subgraph,
                                )?
                            },
                        )
                    } else {
                        None
                    };
                    let coordinate_frames = extend_coordinate_frames(
                        prior_frames,
                        DirectCoordinateFrame {
                            active_subgraph: rescaled_subgraph.clone(),
                            lmb: coordinate_lmb.clone(),
                        },
                    )?;
                    let ctx = UVCtx::new(self.graph, self.settings);
                    // The normalized localization kernel is inert under later Taylor
                    // operations. Transform the complete active CFF, scaling finite
                    // integrated coefficients only when their owners lie in the current
                    // component, and retain the frozen kernel outside the series.
                    let step = Local3DLoopRescaling::FullSubgraph.projection_step(
                        &ctx,
                        current,
                        given,
                        &active_subgraph,
                        Some(&rescaled_subgraph),
                        &canonical_lmb,
                        &coordinate_lmb,
                    );
                    Ok(DirectSector {
                        active_subgraph,
                        coordinate_frames,
                        projection_paths: paths_with_step(&prior_paths, &step),
                        active: -apply_taylor(
                            &ctx,
                            self.localizer.orientation,
                            current,
                            given,
                            Some(rescaled_subgraph.clone()),
                            &coordinate_lmb,
                            active,
                            on_shell.as_ref(),
                        )?,
                        frozen_integrands,
                    })
                },
            )
            .collect::<Result<Vec<_>>>()?;

        let local = Direct3dCts::from_sectors(next)?.map(|atom| {
            Ok(UvMarker::new(self.settings).apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })?;
        crate::debug_tags!(#generation, #uv, #local, #three_d, #profile;
            stage = "direct_3d_construction",
            graph = %self.graph.name,
            current = %current.log_display(),
            given = %given.log_display(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Constructed direct local-3D Taylor sectors"
        );
        Ok(local)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn run_integrated<S: ForestNodeLike, I: ForestNodeLike, M: ForestNodeLike>(
        &mut self,
        finite_counterterms: &[(&Atom, &I)],
        integrated_node: &I,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Direct3dCts> {
        // Each finite effective vertex retains its connected source owner.
        // Equal masses in disjoint integrated prefixes must remain distinct
        // until an enclosing operation contains their respective components.
        let finite_counterterm =
            finite_counterterms
                .iter()
                .fold(Atom::one(), |product, (coefficient, component)| {
                    let mass_scope =
                        function!(*LOCAL_3D_MASS_SCOPE, usize::from(component.lmb_id()) as i64);
                    // Preserve each completed component's vacuum mass before multiplying
                    // coefficients. A disjoint Taylor operation leaves it fixed;
                    // an enclosing operation gives it the same hard weight as
                    // the component's loop momenta. The one-argument mUV head
                    // retains the ordinary mass attributes and no Taylor-dependent
                    // payload.
                    let mut coefficient = coefficient
                        .replace(GS.m_uv_vacuum)
                        .with(function!(GS.m_uv_vacuum, component.subgraph().symbol()));
                    let mut masses = Vec::new();
                    for (pair, _, edge) in self.graph.iter_edges_of(component.subgraph()) {
                        let mass = edge.data.mass_atom();
                        if pair.is_paired() && !mass.is_zero() && !masses.contains(&mass) {
                            // The soft jet holds scalar mass logarithms fixed.
                            // Only algebraic mass factors carry the effective
                            // vertex's engineering weight in the hard-dual chart.
                            coefficient = coefficient
                                .replace(mass.clone())
                                .max_level(0)
                                .with(&mass * &mass_scope);
                            masses.push(mass);
                        }
                    }
                    product * coefficient
                });
        let (active, frozen_integrands) =
            self.localize_integrated(&finite_counterterm, integrated_node)?;
        let active_subgraph = current.reduced_subgraph(given);
        let (canonical_lmb, coordinate_lmb) = coordinate_lmb(
            &UVCtx::new(self.graph, self.settings),
            current,
            given,
            None,
            &[],
            &active_subgraph,
        )?;
        let on_shell = if current.renormalization_scheme() == crate::uv::ApproximationType::OS {
            Some(
                if self.graph.paired_edges(current.subgraph())
                    == self.graph.paired_edges(&self.source_frame.active_subgraph)
                {
                    self.on_shell_source(
                        current,
                        given,
                        &coordinate_lmb,
                        &active,
                        &active_subgraph,
                    )?
                } else if current.subgraph().includes(integrated_node.subgraph())
                    && given.subgraph() == integrated_node.subgraph()
                {
                    if finite_counterterms
                        .iter()
                        .any(|(_, owner)| !integrated_node.subgraph().includes(owner.subgraph()))
                    {
                        return Err(eyre!(
                            "a nested finite source requires its explicit contained owners"
                        ));
                    }
                    self.nested_finite_on_shell_source(
                        current,
                        given,
                        integrated_node,
                        &finite_counterterm,
                        &canonical_lmb,
                        &coordinate_lmb,
                        &active_subgraph,
                    )?
                } else {
                    if integrated_node.subgraph().intersects(current.subgraph())
                        || !given.subgraph().is_empty()
                        || current.dod() != 1
                        || self.graph.n_loops(current.subgraph()) != 1
                        || finite_counterterms.iter().any(|(_, owner)| {
                            owner.subgraph().intersects(current.subgraph())
                                || !integrated_node.subgraph().includes(owner.subgraph())
                        })
                    {
                        return Err(eyre!(
                            "proper direct OS finite replay requires explicit disjoint integrated owners"
                        ));
                    }
                    // The given analytic coefficients and their physical/vacuum
                    // mass scopes remain outside this component's fixed point.
                    // Its normalized localization kernel stays in `frozen_integrands`.
                    let coefficient = self.on_shell_component(current, given, &coordinate_lmb)?;
                    self.attach_on_shell_components(
                        &integrated_node.subgraph().union(current.subgraph()),
                        &coefficient,
                        &finite_counterterm,
                    )?
                },
            )
        } else {
            None
        };
        let ctx = UVCtx::new(self.graph, self.settings);
        let active = -apply_taylor(
            &ctx,
            self.localizer.orientation,
            current,
            given,
            Some(active_subgraph.clone()),
            &coordinate_lmb,
            &active,
            on_shell.as_ref(),
        )?;
        let step = Local3DLoopRescaling::ReducedSubgraph.projection_step(
            &ctx,
            current,
            given,
            &active_subgraph,
            Some(&active_subgraph),
            &canonical_lmb,
            &coordinate_lmb,
        );
        let sector = DirectSector {
            active_subgraph: active_subgraph.clone(),
            projection_paths: paths_with_step(&[], &step),
            coordinate_frames: vec![DirectCoordinateFrame {
                active_subgraph,
                lmb: coordinate_lmb,
            }],
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

#[cfg(test)]
#[path = "source_owner_tests.rs"]
mod source_owner_tests;

#[cfg(test)]
#[path = "nested_source_tests.rs"]
mod nested_source_tests;
