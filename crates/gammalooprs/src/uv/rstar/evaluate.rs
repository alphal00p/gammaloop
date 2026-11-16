//! Assemble the ordinary UV preparation and the UV form of the IR correction.

use std::{
    collections::BTreeMap,
    sync::Arc,
    time::{Duration, Instant},
};

use eyre::{Result, WrapErr, ensure, eyre};
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{SuBitGraph, SubSetLike},
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    parse, symbol,
};
use vakint::{Vakint, VakintSettings};

use crate::{
    debug_tags,
    graph::Graph,
    model::Model,
    utils::GS,
    uv::{UVgenerationSettings, UltravioletGraph, uv_graph::UVE},
};

use super::{
    ChildMass, CountertermCall, IrForestTerm, MassiveUvEvaluator, RStarPlan, UvChildrenRoute,
    UvPreparationTerm,
    boundary::{self, BoundaryMode, BoundaryRequest, Built, ClosedRootNormalization},
    expression::{IntegrationTimings, RStarIntegrand, RootProjector},
    moments::hard_moments,
    physical::PreparedPhysicalUvTerm,
};

impl RStarPlan {
    /// Complete every original proper-UV child and prepare its physical
    /// quotient without reducing or integrating the remaining loop tensors.
    ///
    /// The optional linear root projector acts on each complete quotient
    /// numerator exactly once. Its coefficients must be independent of the
    /// remaining loop momenta. Child Laurent products and their correlations
    /// remain inside the terms, including when an independently grouped bare
    /// amplitude vanishes. Group matching terms with [`super::PhysicalUvGroup`]
    /// and apply the final `-K` to the resulting preparation plus bare poles.
    pub fn prepare_physical_uv_terms(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_projector: Option<&RootProjector<'_>>,
    ) -> Result<Vec<PreparedPhysicalUvTerm>> {
        let evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            None,
            ChildMass::Zero,
        )?;
        let started = Instant::now();
        let settings: Arc<str> = serde_json::to_string(&uv_settings.vakint)?.into();
        let terms = self
            .ordinary_uv
            .iter()
            .filter(|term| !term.counterterms.is_empty())
            .map(|term| {
                PreparedPhysicalUvTerm::new(
                    evaluation.physical_remainder(term)?,
                    model,
                    Arc::clone(&settings),
                    root_projector,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        tracing::info!(
            uv = true, rstar = true, timing = true,
            stage = "rstar_prepare_physical_quotients",
            graph = %graph.name,
            prepared_quotients = terms.len(),
            proper_uv_seconds = evaluation.proper_uv_elapsed.as_secs_f64(),
            preparation_seconds = started.elapsed().as_secs_f64(),
            "Prepared original-forest physical quotients before integration"
        );
        Ok(terms)
    }

    /// Compute the local MS counterterm of a homogeneous massless graph.
    ///
    /// The result retains its open tensor indices and external edge momenta.
    /// Global graph prefactors and external projectors are not included.
    /// All proper UV maps use the existing massive R operation; only the
    /// remaining single-mass integrals use Vakint's p-integral backend. Include
    /// `forcer` in the evaluation methods for that backend.
    pub fn counterterm(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            None,
            ChildMass::Reference,
        )?
        .counterterm()
    }

    /// [`Self::counterterm`] with the proper UV children supplied by the caller
    /// (the polynomials `evaluate_uv_children` would return for this graph and
    /// this plan, keyed by support, with the reference mass
    /// `sqrt(rstar_reference_mass_squared)` on the marked line). Lets a caller
    /// that evaluates many related graphs reuse identical children.
    pub fn counterterm_with_children(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        children: BTreeMap<SuBitGraph, Atom>,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            Some(children),
            ChildMass::Reference,
        )?
        .counterterm()
    }

    /// Compute a root-projected local MS counterterm using the same R* forest.
    ///
    /// The linear projector acts on each complete outer Taylor polynomial,
    /// before its remaining vacuum tensor reduction. It is never applied to
    /// proper UV children or decorated IR quotient counterterms, whose open
    /// tensor boundaries are retained. Its coefficients must be regular in
    /// epsilon and it must not introduce loop momenta. Rational functions of
    /// external momenta are allowed when their denominators are scalar and
    /// independent of every loop momentum. These factors are spectators:
    /// the root Taylor operation has already finished.
    ///
    /// Root normalization factors such as `1/(d-1)` participate in the final
    /// Laurent projection. Scalar marker coefficients may carry several root
    /// projection channels through a single evaluation. Graph-wide prefactors
    /// remain excluded, as for [`Self::counterterm`].
    pub fn projected_counterterm(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_projector: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_projector),
            None,
            ChildMass::Reference,
        )?
        .counterterm()
    }

    /// Local export patch: run only the ordinary-UV preparation terms of
    /// [`Self::projected_counterterm`] that carry no proper UV counterterm
    /// (the empty spinney), through the same root projector, outer Taylor
    /// expansion and single-mass integration. Proper UV children and IR
    /// corrections are not computed. Returns the number of terms integrated;
    /// the result values are discarded (meaningful only with requests skipped
    /// through `PINT6_SKIP_MIN_P_LOOPS`, for request export).
    pub fn projected_bare_terms_only(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_projector: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<usize> {
        self.bare_terms_only_with(graph, model, vakint, uv_settings, Some(root_projector))
    }

    /// As [`Self::projected_counterterm`] for a metric root projector whose
    /// momentum-free metric is applied in the root scope instead: every
    /// root-level remainder is multiplied by [`RStarIntegrand::root_closure`]
    /// together with its child insertions, so the remainder algebra contracts
    /// a closed scalar. IR quotient counterterms and proper UV children keep
    /// their open indices; IR corrections are closed after their quotient `K`.
    /// `root_normalization` then supplies only the momentum-dependent rest of
    /// the projector (normalization, rerouting, channel tags) after Taylor.
    pub fn projected_counterterm_closing_root(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            None,
            ChildMass::Reference,
        )?
        .closing_root()?
        .counterterm()
    }

    /// As [`Self::projected_counterterm_from_pole_part`], closing the root
    /// metric in the root scope (see [`Self::projected_counterterm_closing_root`]).
    /// `bare_poles` keep their open slots; they are closed where the metric
    /// projector applied before, after their own pole projection.
    pub fn projected_counterterm_from_pole_part_closing_root(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        bare_poles: &Atom,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            None,
            ChildMass::Zero,
        )?
        .closing_root()?
        .counterterm_from_pole_part(bare_poles)
    }

    /// Local export patch: as [`Self::projected_bare_terms_only`], closing the
    /// metric root in the root scope (see [`Self::projected_counterterm_closing_root`]).
    pub fn projected_bare_terms_only_closing_root(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<usize> {
        self.bare_terms_only_closing(
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            true,
            None,
        )
    }

    /// [`Self::projected_counterterm_closing_root`] with the closed metric-root
    /// normalisation also given as data (`closed_root`: the callback is
    /// `closed_root.factor · σ` on closed scalars). With
    /// `GAMMALOOP_REMAINDER_BOUNDARY` unset this is
    /// [`Self::projected_counterterm_closing_root`]; otherwise the admitted
    /// root-scope remainders take the remainder-boundary path.
    pub fn projected_counterterm_closing_root_scalar(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
        closed_root: &ClosedRootNormalization,
    ) -> Result<Atom> {
        let mut evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            None,
            ChildMass::Reference,
        )?
        .closing_root()?;
        evaluation.closed_root = Some(closed_root);
        evaluation.counterterm()
    }

    /// [`Self::projected_bare_terms_only_closing_root`] with the closed
    /// metric-root normalisation as data (see
    /// [`Self::projected_counterterm_closing_root_scalar`]).
    pub fn projected_bare_terms_only_closing_root_scalar(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        root_normalization: &dyn Fn(&Atom) -> Result<Atom>,
        closed_root: &ClosedRootNormalization,
    ) -> Result<usize> {
        self.bare_terms_only_closing(
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_normalization),
            true,
            Some(closed_root),
        )
    }

    /// Local export patch: as [`Self::projected_bare_terms_only`], without a
    /// root projector, i.e. the bare terms of [`Self::pole_part`] (the route
    /// of insertion surrogates above four loops).
    pub fn bare_terms_only_unprojected(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<usize> {
        self.bare_terms_only_with(graph, model, vakint, uv_settings, None)
    }

    fn bare_terms_only_with<'a>(
        &'a self,
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        uv_settings: &'a UVgenerationSettings,
        root_projector: Option<&'a RootProjector<'a>>,
    ) -> Result<usize> {
        self.bare_terms_only_closing(
            graph,
            model,
            vakint,
            uv_settings,
            root_projector,
            false,
            None,
        )
    }

    fn bare_terms_only_closing<'a>(
        &'a self,
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        uv_settings: &'a UVgenerationSettings,
        root_projector: Option<&'a RootProjector<'a>>,
        close_root: bool,
        closed_root: Option<&'a ClosedRootNormalization>,
    ) -> Result<usize> {
        let mut evaluation =
            // Bare terms only: no proper UV child is computed.
            Evaluation::new_with(
                self,
                graph,
                model,
                vakint,
                uv_settings,
                root_projector,
                None,
                ChildMass::Reference,
                false,
            )?;
        if close_root {
            evaluation = evaluation.closing_root()?;
        }
        evaluation.closed_root = closed_root;
        let mut count = 0;
        for term in &self.ordinary_uv {
            if !term.counterterms.is_empty() || term.remainder.is_scaleless(graph, self.marked_edge)
            {
                continue;
            }
            super::pint6_dump::set_term(format!(
                "graph={} marked_edge={} term=ordinary_uv[empty spinney] remaining_loops={} taylor_order={}",
                graph.name,
                usize::from(self.marked_edge),
                term.remainder.loops,
                self.taylor_order
            ));
            let _ = evaluation.integrate_ordinary(term)?;
            count += 1;
        }
        evaluation.log_boundary_summary();
        Ok(count)
    }

    /// The local MS counterterm of a logarithmically divergent graph whose
    /// open slots are four Minkowski slots of one representation, through its
    /// three scalar metric channels instead of the open rank-four tensor.
    ///
    /// The counterterm is `sum_j K[P_j . X] T_j` over the metric pairings
    /// `T_j` ([`RStarIntegrand::isotropic_boundary`]) with the dual projectors
    /// `P_j = ((d+2) T_j - sum_i T_i) / (d (d-1) (d+2))`: MS acts on the
    /// coefficients of a d-dimensional basis, which these projections are.
    /// The momentum-free polynomial part of `sum_j tau_j P_j` closes every
    /// root-scope remainder (as [`Self::projected_counterterm_closing_root`]);
    /// its rational normalization joins after the Taylor step, inside this
    /// graph's `K`. Nested children and quotients keep their open slots. The
    /// channel markers `tau_j` are then replaced by `T_j`.
    pub fn isotropic_counterterm(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        pairings: &[Atom; 3],
    ) -> Result<Atom> {
        let d = Atom::var(GS.dim);
        let markers = [
            symbol!("gammalooprs::uv_isotropic_channel_1"; Scalar),
            symbol!("gammalooprs::uv_isotropic_channel_2"; Scalar),
            symbol!("gammalooprs::uv_isotropic_channel_3"; Scalar),
        ];
        // sum_j tau_j ((d+2) T_j - sum_i T_i) = sum_j ((d+2) tau_j - sum_i tau_i) T_j.
        let markers_sum = markers
            .iter()
            .fold(Atom::Zero, |sum, marker| sum + Atom::var(*marker));
        let closure =
            markers
                .iter()
                .zip(pairings)
                .fold(Atom::Zero, |closure, (marker, pairing)| {
                    closure + ((&d + 2) * Atom::var(*marker) - &markers_sum) * pairing
                });
        let normalization = Atom::one() / (&d * (&d - 1) * (&d + 2));
        let normalize = |tensor: &Atom| -> Result<Atom> { Ok(tensor * &normalization) };
        let mut evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(&normalize),
            None,
            ChildMass::Reference,
        )?;
        evaluation.root_closure = Some(closure);
        evaluation.channel = ("isotropic", "1/(d (d-1) (d+2))");
        let channels = evaluation.counterterm()?;
        Ok(markers
            .iter()
            .zip(pairings)
            .fold(channels, |tensor, (marker, pairing)| {
                tensor.replace(Atom::var(*marker)).with(pairing.clone())
            }))
    }

    /// Reconstruct all poles of the bare, infrared-finite massless propagator
    /// graph at its original off-shell momentum, including subdivergence poles.
    ///
    /// From `Delta_U(G) = -K R'_U G`, this returns
    /// `K G = -Delta_U(G) - K sum_{S != empty} Delta_U(S) G/S`.
    /// Proper UV maps use the existing massive R operation; their contracted
    /// remainders are lower-loop p-integrals, with no Taylor expansion or
    /// auxiliary mass. Vakint must support these lower-loop integrals and the
    /// extra epsilon orders required by the insertion poles. Physical IR
    /// divergences in the input are not reconstructed by this identity.
    ///
    /// As for `counterterm`, graph prefactors and external projectors are
    /// excluded; open tensor indices and external momenta are retained.
    pub fn pole_part(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<Atom> {
        let (_, external_routing) = RStarIntegrand::physical_external_routing(
            graph,
            &graph.external_filter::<SuBitGraph>(),
        )?;
        let mut evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            None,
            ChildMass::Reference,
        )?;
        let result = -evaluation.counterterm()? - evaluation.physical_uv_preparation()?;
        evaluation.poles(&result.replace_multiple(&external_routing))
    }

    /// Recover the original local MS counterterm from independently known bare
    /// poles: `Delta_U(G) = -K [K G + sum_{S != empty} Delta_U(S) G/S]`.
    ///
    /// This retains the original graph's proper UV forest. It therefore also
    /// accepts bare poles obtained by integrating propagator insertions or
    /// replacing them with normalized scalar bubbles: those substitutions need
    /// not preserve the original subdivergences or their nested pole maps.
    ///
    /// `bare_poles` must describe this infrared-finite massless two-point graph,
    /// in the supplied Vakint normalization, with its original open tensor
    /// slots and external momentum labels. Exclude graph-wide prefactors and
    /// external projectors, as in [`Self::pole_part`]. Boundary momenta are
    /// identified using the same physical routing as that method. All proper
    /// UV counterterms belong to the original graph; the root R*
    /// counterterm is not evaluated.
    pub fn counterterm_from_pole_part(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        bare_poles: &Atom,
    ) -> Result<Atom> {
        // Only the physical remainders consume the children, at zero mass.
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            None,
            ChildMass::Zero,
        )?
        .counterterm_from_pole_part(bare_poles)
    }

    /// Project the complete physical preparation before its final pole map.
    ///
    /// This is [`Self::counterterm_from_pole_part`] with a linear root projector
    /// whose coefficients are regular in epsilon and contain no loop momenta.
    /// Proper UV children retain their open tensor boundaries. Contractions
    /// and scalar dimensions introduced by the projector are completed before
    /// `K`, including normalization factors such as `1/(d-1)`.
    pub fn projected_counterterm_from_pole_part(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        bare_poles: &Atom,
        root_projector: &dyn Fn(&Atom) -> Result<Atom>,
    ) -> Result<Atom> {
        Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            Some(root_projector),
            None,
            ChildMass::Zero,
        )?
        .counterterm_from_pole_part(bare_poles)
    }
}

#[cfg(test)]
impl RStarPlan {
    /// Tests: the outer Taylor polynomial of every ordinary-UV remainder of a
    /// closing-root evaluation (children included), formed exactly as
    /// [`RStarIntegrand::integrate`] forms it, with the remainder's loop
    /// count, its contracted expression and the idenso certificate of that
    /// contraction ([`RStarIntegrand::contracted_certified`]). Nothing is
    /// projected or integrated.
    pub(super) fn closing_root_outer_taylor_polynomials(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<Vec<(usize, Atom, bool, Atom)>> {
        let evaluation = Evaluation::new(
            self,
            graph,
            model,
            vakint,
            uv_settings,
            None,
            None,
            ChildMass::Reference,
        )?
        .closing_root()?;
        let mut polynomials = Vec::new();
        for term in &self.ordinary_uv {
            if term.remainder.is_scaleless(graph, self.marked_edge) {
                continue;
            }
            let remainder = evaluation.remainder(term)?;
            let (expression, certified) = remainder.contracted_certified()?;
            let polynomial = remainder.outer_taylor(
                expression,
                &evaluation.root_edges,
                usize::try_from(self.taylor_order)?,
                evaluation.reference_mass_squared,
                &evaluation.auxiliary_mass_squared,
                &mut super::expression::IntegrationTimings::new(),
            )?;
            polynomials.push((
                remainder.region.loops,
                expression.clone(),
                *certified,
                polynomial,
            ));
        }
        Ok(polynomials)
    }
}

struct Evaluation<'a> {
    plan: &'a RStarPlan,
    graph: &'a Graph,
    model: &'a Model,
    vakint: &'a Vakint,
    uv_settings: &'a UVgenerationSettings,
    settings: VakintSettings,
    epsilon: Symbol,
    reference_mass_squared: Symbol,
    auxiliary_mass_squared: Atom,
    children: BTreeMap<SuBitGraph, Atom>,
    /// Carriers scaled by the outer Taylor operation: the graph's external
    /// edges, or the full crown of a proper support root.
    root_edges: Vec<EdgeIndex>,
    uv: MassiveUvEvaluator<'a>,
    quotient_counterterms: BTreeMap<(CountertermCall, Atom), Atom>,
    integration_id: usize,
    root_projector: Option<&'a RootProjector<'a>>,
    /// The root metric closure of a closing evaluation (root scope only).
    root_closure: Option<Atom>,
    /// The root projection channel and its normalisation, recorded in the
    /// closure dump headers (`PINT6_CLOSURE_DUMP`).
    channel: (&'static str, &'static str),
    /// The closed metric-root normalisation as data, given by the `_scalar`
    /// entry points; used only when `GAMMALOOP_REMAINDER_BOUNDARY` is set.
    closed_root: Option<&'a ClosedRootNormalization>,
    /// Time on the remainder-boundary path (inert when the mode is off).
    boundary_elapsed: Duration,
    /// Remainder-boundary counts of this evaluation (`BOUNDARY_SUMMARY`).
    summary: boundary::BoundarySummary,
    started: Instant,
    proper_uv_elapsed: Duration,
    quotient_elapsed: Duration,
    integration_elapsed: Duration,
}

impl<'a> Evaluation<'a> {
    fn new(
        plan: &'a RStarPlan,
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        uv_settings: &'a UVgenerationSettings,
        root_projector: Option<&'a RootProjector<'a>>,
        precomputed_children: Option<BTreeMap<SuBitGraph, Atom>>,
        child_mass: ChildMass,
    ) -> Result<Self> {
        Self::new_with(
            plan,
            graph,
            model,
            vakint,
            uv_settings,
            root_projector,
            precomputed_children,
            child_mass,
            true,
        )
    }

    fn new_with(
        plan: &'a RStarPlan,
        graph: &'a Graph,
        model: &'a Model,
        vakint: &'a Vakint,
        uv_settings: &'a UVgenerationSettings,
        root_projector: Option<&'a RootProjector<'a>>,
        precomputed_children: Option<BTreeMap<SuBitGraph, Atom>>,
        child_mass: ChildMass,
        compute_children: bool,
    ) -> Result<Self> {
        let started = Instant::now();
        for (_, _, edge) in graph.iter_edges_of(&plan.support) {
            ensure!(
                edge.data.mass_atom(model).is_zero(),
                "R* currently requires a massless input graph"
            );
        }
        let settings = uv_settings.vakint.true_settings();
        let epsilon = parse!(&settings.epsilon_symbol)
            .get_symbol()
            .ok_or_else(|| eyre!("Vakint epsilon must be a symbol"))?;
        let reference_mass_squared = *super::REFERENCE_MASS_SQUARED;
        let reference_mass = Atom::var(reference_mass_squared).sqrt();
        let children_started = Instant::now();
        let children = if let Some(children) = precomputed_children {
            children
        } else if plan.counterterms.is_empty() || !compute_children {
            BTreeMap::new()
        } else {
            plan.evaluate_uv_children_for(
                graph,
                model,
                &reference_mass,
                vakint,
                uv_settings,
                child_mass,
            )
            .wrap_err(match plan.uv_children {
                UvChildrenRoute::Massive => {
                    "while computing the proper UV counterterms with massive R"
                }
                UvChildrenRoute::Rstar { .. } => {
                    "while computing the proper UV counterterms with recursive R*"
                }
            })?
        };
        let root_edges = match &plan.root_boundary {
            None => {
                let root_external = graph.external_filter::<SuBitGraph>();
                graph
                    .iter_edges_of(&root_external)
                    .map(|(_, edge, _)| edge)
                    .collect::<Vec<_>>()
            }
            Some(boundary) => graph
                .iter_edges_of(boundary)
                .map(|(_, edge, _)| edge)
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect(),
        };
        let proper_uv_elapsed = children_started.elapsed();
        let mut uv = MassiveUvEvaluator::new(graph, model, vakint, uv_settings);
        if let Some(cache) = &plan.counterterm_cache {
            uv = uv.with_counterterm_cache(cache.clone());
        }
        Ok(Self {
            plan,
            graph,
            model,
            vakint,
            uv_settings,
            settings,
            epsilon,
            reference_mass_squared,
            auxiliary_mass_squared: Atom::var(symbol!("rstar_mass_squared"; Scalar)),
            children,
            root_edges,
            uv,
            quotient_counterterms: BTreeMap::new(),
            integration_id: 0,
            root_projector,
            root_closure: None,
            channel: if root_projector.is_some() {
                (
                    "projection open-root",
                    "root projector after the outer Taylor step",
                )
            } else {
                ("none", "1")
            },
            closed_root: None,
            boundary_elapsed: Duration::ZERO,
            summary: Default::default(),
            started,
            proper_uv_elapsed,
            quotient_elapsed: Duration::ZERO,
            integration_elapsed: Duration::ZERO,
        })
    }

    /// Close the root metric in every root-scope remainder of this evaluation.
    fn closing_root(mut self) -> Result<Self> {
        self.root_closure = Some(RStarIntegrand::root_closure(
            self.graph,
            &self.plan.support,
        )?);
        self.channel = ("projection metric", "1/(i p^2 d_A (d-1))");
        Ok(self)
    }

    fn counterterm(&mut self) -> Result<Atom> {
        if self.graph.compute_dod(&self.plan.support) < 0 {
            self.log_timing();
            return Ok(Atom::Zero);
        }
        let mut prepared = Atom::Zero;
        for term in &self.plan.ordinary_uv {
            if term
                .remainder
                .is_scaleless(self.graph, self.plan.marked_edge)
            {
                continue;
            }
            prepared += self.integrate_ordinary(term)?;
        }
        for correction in &self.plan.ir {
            for term in &correction.preparation {
                if term
                    .remainder
                    .is_scaleless(self.graph, self.plan.marked_edge)
                {
                    continue;
                }
                // The hard part feeds quotient counterterms (a `K` each): it
                // keeps its open root slots, closed after the quotient below.
                let hard = self.open_remainder(term)?;
                let polynomial = hard
                    .taylor(usize::try_from(correction.taylor_order)?)
                    .wrap_err("while Taylor-expanding the prepared hard component")?;
                let covariant = hard.project_hard(&polynomial)?;
                let moments = hard_moments(&covariant, hard.soft_edges()).wrap_err_with(|| {
                    format!(
                        "while separating hard moments: hard={}, remaining={}, soft={:?}, covariant={covariant}",
                        correction.hard.string_label(),
                        term.remainder.remaining.string_label(),
                        hard.soft_edges(),
                    )
                })?;
                for forest in &correction.forest {
                    let mut expression = Atom::Zero;
                    for (coefficient, insertion) in &moments {
                        expression += coefficient * self.ir_forest(forest, insertion)?;
                    }
                    if !expression.is_zero() {
                        if let Some(closure) = &self.root_closure {
                            expression *= closure;
                        }
                        prepared += self.integrate(&hard, &expression)?;
                    }
                }
            }
        }
        let result = -self.poles(&prepared)?;
        self.log_timing();
        Ok(result)
    }

    /// The remainder-boundary summary of this evaluation, when the mode is
    /// not off (`BOUNDARY_SUMMARY`, DESIGN 8).
    fn log_boundary_summary(&self) {
        if !matches!(
            boundary::remainder_boundary_mode(),
            Ok(BoundaryMode::On | BoundaryMode::Check)
        ) {
            return;
        }
        let root = if self.root_closure.is_some() && self.root_projector.is_some() {
            if self.plan.root_boundary.is_some() {
                "child"
            } else {
                "closed"
            }
        } else if self.plan.root_boundary.is_some() || self.root_projector.is_none() {
            "child"
        } else {
            "open"
        };
        self.summary.log(&self.graph.name, root);
    }

    fn log_timing(&self) {
        self.log_boundary_summary();
        let elapsed = self.started.elapsed();
        let preparation = elapsed.saturating_sub(
            self.proper_uv_elapsed + self.quotient_elapsed + self.integration_elapsed,
        );
        tracing::info!(
            uv = true,
            rstar = true,
            timing = true,
            stage = "rstar_counterterm_phase_timing",
            graph = %self.graph.name,
            root_projected = self.root_projector.is_some(),
            total_ms = elapsed.as_millis() as u64,
            proper_uv_ms = self.proper_uv_elapsed.as_millis() as u64,
            ir_quotient_ms = self.quotient_elapsed.as_millis() as u64,
            outer_integration_ms = self.integration_elapsed.as_millis() as u64,
            preparation_ms = preparation.as_millis() as u64,
            quotient_calls = self.quotient_counterterms.len(),
            outer_integrations = self.integration_id,
            "R* counterterm phase timings"
        );
    }

    fn poles(&self, expression: &Atom) -> Result<Atom> {
        Ok(expression.series(self.epsilon, 0, -1)?.to_atom().expand())
    }

    fn counterterm_from_pole_part(&mut self, bare_poles: &Atom) -> Result<Atom> {
        let (_, external_routing) = RStarIntegrand::physical_external_routing(
            self.graph,
            &self.graph.external_filter::<SuBitGraph>(),
        )?;
        // Closed bare poles stand where the metric projector contracted them.
        let bare_poles = match &self.root_closure {
            Some(closure) => bare_poles * closure,
            None => bare_poles.clone(),
        };
        let mut prepared = bare_poles + self.physical_uv_preparation()?;
        if let Some(projector) = self.root_projector {
            prepared = projector(&prepared)?;
        }
        // The independent bare-pole provider and the root projector can both
        // introduce tensor contractions with evanescent scalar coefficients.
        // Complete them before K, just as for an ordinary integrated UV node.
        let prepared = RStarIntegrand::finish_tensor(
            prepared.replace_multiple(&external_routing),
            self.integration_id + 1,
        )?;
        let result = -self.poles(&prepared)?;
        self.log_timing();
        Ok(result)
    }

    /// The physical, nonempty proper-UV sum, before the root pole projection.
    /// Child counterterms remain in their original scopes and multiply their
    /// quotients before the scalar integral's Laurent depth is determined.
    fn physical_uv_preparation(&mut self) -> Result<Atom> {
        let mut result = Atom::Zero;
        for term in self
            .plan
            .ordinary_uv
            .iter()
            .filter(|term| !term.counterterms.is_empty())
        {
            let remainder = self.physical_remainder(term)?;
            self.integration_id += 1;
            let mut settings = self.settings.clone();
            settings.number_of_terms_in_epsilon_expansion = i64::try_from(term.remainder.loops)?;
            let started = Instant::now();
            result += remainder.integrate_physical(
                self.vakint,
                self.uv_settings,
                &settings,
                self.integration_id,
            )?;
            self.integration_elapsed += started.elapsed();
        }
        Ok(result)
    }

    fn physical_remainder(&self, term: &UvPreparationTerm) -> Result<RStarIntegrand> {
        let insertions = self
            .insertions(term)?
            .into_iter()
            .map(|insertion| insertion.replace(self.reference_mass_squared).with(0))
            .collect::<Vec<_>>();
        RStarIntegrand::new_with_closure(
            self.graph,
            self.model,
            &term.remainder,
            &insertions,
            self.plan.marked_edge,
            Atom::Zero,
            self.root_closure.as_ref(),
        )
    }

    fn insertions(&self, term: &UvPreparationTerm) -> Result<Vec<Atom>> {
        term.counterterms
            .iter()
            .map(|child| {
                self.children
                    .get(&child.argument.support)
                    .cloned()
                    .ok_or_else(|| eyre!("missing a completed proper UV counterterm"))
            })
            .collect()
    }

    /// A root-scope remainder, closed when this evaluation closes the root.
    fn remainder(&self, term: &UvPreparationTerm) -> Result<RStarIntegrand> {
        self.remainder_with(term, self.root_closure.as_ref())
    }

    /// A remainder that keeps its open root slots (an IR hard part).
    fn open_remainder(&self, term: &UvPreparationTerm) -> Result<RStarIntegrand> {
        self.remainder_with(term, None)
    }

    fn remainder_with(
        &self,
        term: &UvPreparationTerm,
        closure: Option<&Atom>,
    ) -> Result<RStarIntegrand> {
        let insertions = self.insertions(term)?;
        RStarIntegrand::new_with_closure(
            self.graph,
            self.model,
            &term.remainder,
            &insertions,
            self.plan.marked_edge,
            Atom::var(self.reference_mass_squared),
            closure,
        )
    }

    fn quotient(&mut self, call: &CountertermCall, insertion: &Atom) -> Result<Atom> {
        let key = (call.clone(), insertion.clone());
        if let Some(value) = self.quotient_counterterms.get(&key) {
            return Ok(value.clone());
        }
        let started = Instant::now();
        let result = self
            .uv
            .counterterm(call, std::slice::from_ref(insertion))
            .wrap_err_with(|| {
                format!(
                    "while applying massive R to the decorated quotient {}",
                    call.argument.remaining.string_label()
                )
            })?;
        self.quotient_elapsed += started.elapsed();
        self.quotient_counterterms.insert(key, result.clone());
        Ok(result)
    }

    fn ir_forest(&mut self, forest: &IrForestTerm, insertion: &Atom) -> Result<Atom> {
        let mut rows = vec![(Atom::num(forest.sign as i64), insertion.clone())];
        for (index, call) in forest.steps.iter().enumerate() {
            let mut next = Vec::new();
            for (coefficient, polynomial) in rows {
                let counterterm = self.quotient(call, &polynomial)?;
                if index + 1 == forest.steps.len() {
                    next.push((coefficient, counterterm));
                    continue;
                }
                // The c-operation product is coefficientwise in epsilon.
                // Retain the complete tensor polynomial inside the next UV
                // call: its metric contractions can depend on the dimension.
                for (power, tensor_polynomial) in counterterm.series(self.epsilon, 0, -1)?.terms() {
                    if !tensor_polynomial.is_zero() {
                        next.push((
                            &coefficient * Atom::var(self.epsilon).pow(power),
                            tensor_polynomial.clone(),
                        ));
                    }
                }
            }
            rows = next;
        }
        Ok(rows
            .into_iter()
            .fold(Atom::Zero, |sum, (coefficient, polynomial)| {
                sum + coefficient * polynomial
            }))
    }

    /// Remainder-boundary admission (DESIGN 10a): every reason is silent and
    /// counted, never an error, except the vakint gates under STRICT.
    fn boundary_admits(&self, remainder: &RStarIntegrand) -> std::result::Result<(), String> {
        if self.closed_root.is_none() {
            return Err("no_closed_root".into());
        }
        if self.root_closure.is_none() || self.root_projector.is_none() {
            return Err("open_root".into());
        }
        if self.plan.root_boundary.is_some() {
            return Err("proper_support_root".into());
        }
        if remainder.region.loops < 2 {
            return Err("one_loop_tree_cut".into());
        }
        if !RStarIntegrand::plain_remainder_contraction() {
            return Err("edge_contraction".into());
        }
        if [
            "direct_lowering",
            "tensor_factors",
            "topology_groups",
            "term_rewrites",
            "factor_maps",
            "dot_notation",
        ]
        .iter()
        .any(|name| crate::uv::approx::streamed::legacy(name))
        {
            return Err("algebra_legacy".into());
        }
        if self.settings.verify_numerator_identification {
            return Err("strict_split_branch".into());
        }
        boundary::vakint_gates()
            .and_then(|()| boundary::settings_gate(&self.settings))
            .map_err(|reason| format!("vakint_gates:{reason}"))
    }

    /// One root-scope ordinary remainder (the ordinary loop and the census
    /// bare-term loop): today's [`Self::integrate`] with the mode off; the
    /// remainder-boundary path for an admitted remainder otherwise (DESIGN
    /// 6a, 7).
    fn integrate_ordinary(&mut self, term: &UvPreparationTerm) -> Result<Atom> {
        let mode = boundary::remainder_boundary_mode()?;
        let remainder = self.remainder(term)?;
        if mode == BoundaryMode::Off {
            return self.integrate(&remainder, remainder.expression()?);
        }
        self.summary.remainders += 1;
        let strict = boundary::boundary_strict();
        if let Err(reason) = self.boundary_admits(&remainder) {
            if strict && reason.starts_with("vakint_gates") {
                return Err(eyre!(
                    "BOUNDARY strict: the vakint gates reject an admissible remainder of {} ({reason})",
                    self.graph.name
                ));
            }
            *self.summary.not_admitted.entry(reason).or_default() += 1;
            return self.integrate(&remainder, remainder.expression()?);
        }
        let closed_root = self
            .closed_root
            .expect("admission requires the closed-root normalisation");
        let id = self.integration_id + 1;
        let started = Instant::now();
        let mut settings = self.settings.clone();
        settings.number_of_terms_in_epsilon_expansion = i64::try_from(remainder.region.loops)?;
        let degree = usize::try_from(self.plan.taylor_order)?;
        let built = boundary::build_boundary(
            &remainder,
            &BoundaryRequest {
                graph_name: &self.graph.name,
                marked_edge: self.plan.marked_edge,
                integration_id: id,
                root_edges: &self.root_edges,
                degree,
                reference_mass_squared: self.reference_mass_squared,
                auxiliary_mass_squared: &self.auxiliary_mass_squared,
                vakint: self.vakint,
                settings: &settings,
                additional_normalization: &self.uv_settings.vakint.additional_normalization,
                epsilon: self.epsilon,
                closed_root,
            },
        )
        .wrap_err_with(|| {
            format!(
                "while building the remainder boundary {id} of {}",
                self.graph.name
            )
        })?;
        let large_limit = boundary::fallback_max_summands();
        let built = match built {
            Built::Declined(decline, stats) => {
                let large = boundary::large_decline(&stats, remainder.region.loops);
                eprintln!(
                    "BOUNDARY_FALLBACK graph={} id={id} loops={} stage={} reason={} p_summands={} large={large}",
                    self.graph.name,
                    remainder.region.loops,
                    decline.stage,
                    decline.reason,
                    stats.p_summands
                );
                self.summary.declined(&decline);
                self.summary.contracted(&stats);
                self.boundary_elapsed += started.elapsed();
                self.summary.seconds_boundary += started.elapsed().as_secs_f64();
                if large {
                    self.summary.errors += 1;
                    return Err(eyre!(
                        "BOUNDARY_DECLINED_LARGE graph={} id={id} loops={}: {}:{} ({}; no fallback to today's path)",
                        self.graph.name,
                        remainder.region.loops,
                        decline.stage,
                        decline.reason,
                        if stats.contracted {
                            format!("{} summands of P, limit {large_limit}", stats.p_summands)
                        } else {
                            format!(
                                "before the contraction, at least {} loops or a census dump",
                                boundary::LARGE_PRE_CONTRACTION_LOOPS
                            )
                        }
                    ));
                }
                if strict {
                    self.summary.errors += 1;
                    return Err(eyre!(
                        "BOUNDARY strict: the remainder boundary {id} of {} declined at {}: {}",
                        self.graph.name,
                        decline.stage,
                        decline.reason
                    ));
                }
                return self.integrate(&remainder, remainder.expression()?);
            }
            built => built,
        };
        match mode {
            BoundaryMode::Check => self.check_boundary(&remainder, built, &settings, id, started),
            _ => self.evaluate_boundary(&remainder, built, &settings, id, started, large_limit),
        }
    }

    /// Mode on: the value of a built boundary, finished as today's outer
    /// integration (one integration id, kind `single_mass`).
    fn evaluate_boundary(
        &mut self,
        remainder: &RStarIntegrand,
        built: Built,
        settings: &VakintSettings,
        id: usize,
        started: Instant,
        large_limit: usize,
    ) -> Result<Atom> {
        let loops = remainder.region.loops;
        let (value, mut stats, built) = match built {
            Built::Zero(stats) => {
                self.summary.zero += 1;
                (Atom::Zero, stats, None)
            }
            Built::Boundary(built, stats) => (Atom::Zero, stats, Some(built)),
            Built::Declined(..) => unreachable!("declines are handled by the caller"),
        };
        self.summary.contracted(&stats);
        let mut value = value;
        let evaluate_started = Instant::now();
        if let Some(built) = &built {
            self.summary.boundary += 1;
            let p_loops = built.loops.saturating_sub(1);
            let skipped =
                super::pint6_dump::settings().is_some_and(|dump| p_loops >= dump.skip_min_p_loops);
            if !skipped {
                let options = vakint::p_integral::BoundaryOptions {
                    strict: boundary::boundary_strict() || stats.p_summands > large_limit,
                };
                let (evaluated, reads, hits) =
                    boundary::evaluate_boundary(&built.boundary, self.vakint, settings, options)
                        .wrap_err_with(|| {
                            format!(
                                "while evaluating the remainder boundary {id} of {}",
                                self.graph.name
                            )
                        })?;
                value = evaluated;
                stats.block_reads = reads;
                stats.block_hits = hits;
            }
        }
        let evaluate_seconds = evaluate_started.elapsed().as_secs_f64();
        stats.seconds.push(("evaluate", evaluate_seconds));
        stats.rss_gb[3] = boundary::rss_gb();
        let finish_started = Instant::now();
        let result = RStarIntegrand::finish_tensor(value, id)?;
        let finish_seconds = finish_started.elapsed().as_secs_f64();
        stats.seconds.push(("finish", finish_seconds));
        self.integration_id = id;
        let stage = |name: &str| {
            stats
                .seconds
                .iter()
                .filter(|(stage, _)| *stage == name)
                .map(|(_, seconds)| seconds)
                .sum::<f64>()
        };
        let taylor = stage("den") + stage("collect");
        let conversion = [
            "frames", "cut", "routing", "contract", "project", "assemble", "verify",
        ]
        .iter()
        .map(|name| stage(name))
        .sum::<f64>();
        IntegrationTimings::boundary(
            started,
            taylor,
            conversion,
            evaluate_seconds,
            finish_seconds,
        )
        .log(self.graph, id, loops, "single_mass");
        boundary::log_stats(&self.graph.name, id, &stats);
        boundary::record_stats(&self.graph.name, id, &stats);
        let elapsed = started.elapsed();
        self.integration_elapsed += elapsed;
        self.boundary_elapsed += elapsed;
        self.summary.seconds_boundary += elapsed.as_secs_f64();
        Ok(result)
    }

    /// Mode check (DESIGN 7): the boundary's value, then today's path under
    /// capture (its closures, value and root callback), compared per vacuum
    /// integral; one `BOUNDARY_CHECK` line. Returns today's value; any
    /// MISMATCH is an error.
    fn check_boundary(
        &mut self,
        remainder: &RStarIntegrand,
        built: Built,
        settings: &VakintSettings,
        id: usize,
        started: Instant,
    ) -> Result<Atom> {
        let (built, mut stats) = match built {
            Built::Zero(stats) => {
                self.summary.zero += 1;
                (None, stats)
            }
            Built::Boundary(built, stats) => {
                self.summary.boundary += 1;
                (Some(built), stats)
            }
            Built::Declined(..) => unreachable!("declines are handled by the caller"),
        };
        self.summary.contracted(&stats);
        let closed_root = self
            .closed_root
            .expect("admission requires the closed-root normalisation");
        let p_loops = remainder.region.loops.saturating_sub(1);
        let skipped =
            super::pint6_dump::settings().is_some_and(|dump| p_loops >= dump.skip_min_p_loops);
        let evaluate_started = Instant::now();
        let new_value = match &built {
            Some(_) if skipped => None,
            Some(built) => {
                let options = vakint::p_integral::BoundaryOptions {
                    strict: boundary::boundary_strict()
                        || stats.p_summands > boundary::fallback_max_summands(),
                };
                let (value, reads, hits) =
                    boundary::evaluate_boundary(&built.boundary, self.vakint, settings, options)
                        .wrap_err_with(|| {
                            format!(
                                "while evaluating the remainder boundary {id} of {}",
                                self.graph.name
                            )
                        })?;
                stats.block_reads = reads;
                stats.block_hits = hits;
                Some(value)
            }
            None => Some(Atom::Zero),
        };
        stats
            .seconds
            .push(("evaluate", evaluate_started.elapsed().as_secs_f64()));
        stats.rss_gb[3] = boundary::rss_gb();
        let seconds_new = started.elapsed().as_secs_f64();
        self.summary.seconds_boundary += seconds_new;
        self.boundary_elapsed += started.elapsed();
        // Today's path, captured.
        let old_started = Instant::now();
        let ((today, callback), capture) = super::integrate::capturing_closures(|| {
            boundary::with_callback_check(closed_root, || -> Result<Atom> {
                let expression = remainder.expression()?;
                self.integrate(remainder, expression)
            })
        });
        let today = today?;
        let seconds_old = old_started.elapsed().as_secs_f64();
        let check_started = Instant::now();
        let old_value = if capture.closures.iter().all(|closure| closure.evaluated) {
            capture.value.as_ref()
        } else {
            None
        };
        let label = format!("graph={} id={id}", self.graph.name);
        let report = boundary::compare_pairs(
            built.as_ref(),
            &capture.closures,
            new_value.as_ref(),
            old_value,
            boundary::check_pair_max_bytes(),
            &label,
        );
        let seconds_check = check_started.elapsed().as_secs_f64();
        boundary::log_stats(&self.graph.name, id, &stats);
        boundary::record_stats(&self.graph.name, id, &stats);
        eprintln!(
            "BOUNDARY_CHECK graph={} id={id} loops={} patterns={} zero_patterns={} classes={} old_closures={} pairs={} input={} value={} callback={callback} seconds_new={seconds_new:.3} seconds_old={seconds_old:.3} seconds_check={seconds_check:.3}",
            self.graph.name,
            remainder.region.loops,
            stats.patterns,
            stats.zero_patterns,
            stats.classes,
            report.old_closures,
            report.pairs,
            report.input,
            report.value,
        );
        for line in &report.failing {
            eprintln!("{line}");
        }
        self.summary.checks += 1;
        if report.input.starts_with("skipped") {
            self.summary.checks_skipped += 1;
        }
        if report.mismatch(callback) {
            self.summary.errors += 1;
            return Err(eyre!(
                "BOUNDARY_CHECK mismatch graph={} id={id}: input={} value={} callback={callback}",
                self.graph.name,
                report.input,
                report.value
            ));
        }
        Ok(today)
    }

    fn integrate(&mut self, remainder: &RStarIntegrand, expression: &Atom) -> Result<Atom> {
        self.integration_id += 1;
        // All insertion poles and hard projector coefficients are already in
        // this expression. Only summation and the root -K follow, so request
        // poles of the complete product. The integration helper increases the
        // scalar integral's depth by the coefficient's negative valuation.
        let mut settings = self.settings.clone();
        settings.number_of_terms_in_epsilon_expansion = i64::try_from(remainder.region.loops)?;
        debug_tags!(#uv, #integrated, #rstar;
            integration_id = self.integration_id,
            loops = remainder.region.loops,
            log.expression = expression,
            "Integrating a compensated R* remainder"
        );
        let started = Instant::now();
        let channel = super::pint6_dump::ClosureChannel {
            channel: self.channel.0,
            taylor_order: i64::from(self.plan.taylor_order),
            normalisation: self.channel.1,
        };
        let result = super::pint6_dump::with_channel(channel, || {
            remainder.integrate(
                expression,
                self.graph,
                &self.root_edges,
                usize::try_from(self.plan.taylor_order)?,
                self.reference_mass_squared,
                &self.auxiliary_mass_squared,
                self.vakint,
                self.uv_settings,
                &settings,
                self.integration_id,
                self.root_projector,
            )
        });
        self.integration_elapsed += started.elapsed();
        self.summary.seconds_today += started.elapsed().as_secs_f64();
        result
    }
}
