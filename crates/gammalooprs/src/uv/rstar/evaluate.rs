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
    uv::{UVgenerationSettings, UltravioletGraph, uv_graph::UVE},
};

use super::{
    ChildMass, CountertermCall, IrForestTerm, MassiveUvEvaluator, RStarPlan, UvChildrenRoute,
    UvPreparationTerm,
    expression::{RStarIntegrand, RootProjector},
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
        )?
        .closing_root()?
        .counterterm_from_pole_part(bare_poles)
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
        } else if plan.counterterms.is_empty() {
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
            let remainder = self.remainder(term)?;
            prepared += self.integrate(&remainder, remainder.expression()?)?;
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

    fn log_timing(&self) {
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
        let result = remainder.integrate(
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
        );
        self.integration_elapsed += started.elapsed();
        result
    }
}
