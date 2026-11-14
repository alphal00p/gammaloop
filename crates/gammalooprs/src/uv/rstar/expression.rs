//! Explicit hard Taylor coefficients with the existing UV tensor conventions.

use std::{
    borrow::Cow,
    collections::BTreeSet,
    sync::{LazyLock, OnceLock},
    time::Instant,
};

use eyre::{Result, WrapErr, ensure, eyre};
use feynkit_tensor::TensorReducer;
use idenso::{
    CookMode, CookSettings,
    color::ColorSimplifySettings,
    dirac::GammaSimplifySettings,
    tensor::{AlgebraSettings, ContractSettings, SymbolicTensor},
};
use linnet::half_edge::{
    involution::{EdgeIndex, HedgePair},
    subgraph::{InternalSubGraph, SuBitGraph, SubSetOps},
};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        partial::PartialStructureExt,
        representation::{Minkowski, RepName},
        slot::{SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::Replacement,
    symbol,
};
use vakint::{Vakint, VakintSettings, stage_timing::StageSpan, symbols::S};

use crate::{
    debug_tags,
    graph::{Autogen, Graph, LMBext, LoopMomentumBasis},
    model::Model,
    utils::{GS, W_},
    uv::{
        UVgenerationSettings, UltravioletGraph,
        approx::{
            integrated::{Integrated, to_vakint_integrand},
            streamed,
        },
        uv_graph::UVE,
    },
};

use super::{GraphRegion, TaylorExpansion, integrate_single_mass, project::project_hard_vacuum};

pub(super) type RootProjector<'a> = dyn Fn(&Atom) -> Result<Atom> + 'a;

// The shared transport knows explicit edge momenta and these FeynKit basis
// coordinates. Other live momentum conventions must keep their own occurrence.
static UNSUPPORTED_PHYSICAL_MOMENTA: LazyLock<Vec<Symbol>> = LazyLock::new(|| {
    let supported = [
        GS.emr_mom,
        feynkit_graph::symbols::loop_momentum(),
        feynkit_graph::symbols::external_momentum(),
    ];
    [GS.loop_mom, GS.external_mom, S.k, S.p]
        .into_iter()
        .filter(|head| !supported.contains(head))
        .collect()
});

static MEMORY_RELEASE: std::sync::OnceLock<fn()> = std::sync::OnceLock::new();

/// Set by the host binary to return freed allocator memory to the kernel at
/// once; used before the hard Taylor fallback (`simplify_hard_taylor`).
pub fn set_memory_release(release: fn()) {
    let _ = MEMORY_RELEASE.set(release);
}

/// Exclusive wall-clock stages within one outer integral. Backend time may
/// contain its own finer timing events; those must not be added to this total.
struct IntegrationTimings {
    started: Instant,
    last: Instant,
    routing: f64,
    taylor: f64,
    materialization: f64,
    tensor_projection: f64,
    algebra: f64,
    conversion: f64,
    backend: f64,
    finish: f64,
}

impl IntegrationTimings {
    fn new() -> Self {
        let started = Instant::now();
        Self {
            started,
            last: started,
            routing: 0.0,
            taylor: 0.0,
            materialization: 0.0,
            tensor_projection: 0.0,
            algebra: 0.0,
            conversion: 0.0,
            backend: 0.0,
            finish: 0.0,
        }
    }

    fn lap(&mut self) -> f64 {
        let now = Instant::now();
        let elapsed = now.duration_since(self.last).as_secs_f64();
        self.last = now;
        elapsed
    }

    fn log(&self, graph: &Graph, integration_id: usize, loops: usize, kind: &str) {
        tracing::info!(
            uv = true,
            rstar = true,
            timing = true,
            stage = "rstar_outer_integration_detail",
            graph = %graph.name,
            integration_id,
            loops,
            kind,
            total_seconds = self.started.elapsed().as_secs_f64(),
            routing_seconds = self.routing,
            taylor_seconds = self.taylor,
            materialization_seconds = self.materialization,
            tensor_projection_seconds = self.tensor_projection,
            algebra_seconds = self.algebra,
            conversion_seconds = self.conversion,
            backend_seconds = self.backend,
            finish_seconds = self.finish,
            "Outer integration stage timings"
        );
    }
}

/// `GAMMALOOP_HARD_TAYLOR_CHECK=1` also runs the contracted hard Taylor path on
/// every hard component and reports `HARD_TAYLOR_CHECK ... equal=...` on stderr.
fn hard_taylor_check() -> bool {
    static CHECK: LazyLock<bool> =
        LazyLock::new(|| std::env::var_os("GAMMALOOP_HARD_TAYLOR_CHECK").is_some());
    *CHECK
}

/// The original tensor slots and loop carriers survive all graph contractions.
/// Only propagator endpoints are identified in the converter's graph.
pub(super) struct RStarIntegrand {
    pub region: GraphRegion,
    /// The routed product of the remaining numerator, insertions and
    /// denominators, before any index contraction. A hard component is only
    /// ever Taylor projected, and its truncation keeps a tiny fraction of the
    /// contracted terms: it is contracted after the Taylor step.
    routed: Atom,
    /// [`Self::routed`] contracted to dots, built on first use.
    contracted: OnceLock<Atom>,
    graph: Graph,
    basis: LoopMomentumBasis,
    soft_edges: Vec<EdgeIndex>,
}

impl RStarIntegrand {
    pub fn new(
        graph: &Graph,
        model: &Model,
        region: &GraphRegion,
        insertions: &[Atom],
        marked: EdgeIndex,
        reference_mass_squared: Atom,
    ) -> Result<Self> {
        let _stage = StageSpan::new("gamma.remainder.construct");
        ensure!(
            insertions.len() == region.insertions.len(),
            "each R* contraction requires an insertion polynomial"
        );
        let mut contracted: SuBitGraph = graph.empty_subgraph();
        for slot in &region.insertions {
            contracted.union_with(&slot.support);
        }
        ensure!(
            region.support.subtract(&contracted) == region.remaining,
            "R* remainder does not match its contracted supports"
        );
        let numerator_stage = StageSpan::new("gamma.remainder.numerator");
        let mut expression = graph
            .numerator(&region.remaining, &contracted)
            .to_d_dim(GS.dim)
            .get_single_atom()?;
        drop(numerator_stage);
        let routing_stage = StageSpan::new("gamma.remainder.routing_and_denominators");
        for insertion in insertions {
            expression *= insertion;
        }
        // Diagnostics: the remainder before routing (edge momenta intact),
        // numerator times insertions; its routed form is `hard_routed`.
        streamed::dump("hard_source", expression.as_view());
        let contracted = InternalSubGraph::try_new(contracted, &graph.underlying)
            .ok_or_else(|| eyre!("an R* insertion contains an unpaired half-edge"))?;
        let basis = if region.loops == 0 {
            graph.empty_lmb()
        } else {
            graph
                .shrunken_sub_lmb(&region.support, &contracted, region.boundary.clone())
                .wrap_err("could not route an R* prepared component")?
        };
        let soft_edges = graph
            .iter_edges_of(&region.boundary)
            .map(|(_, edge, _)| edge)
            .collect();
        let experiment_input = (super::edge_experiment::enabled() && region.loops != 0)
            .then(|| expression.clone());
        if region.loops != 0 {
            // Boundary Q labels are the independent source variables of the
            // local polynomial. A subgraph basis supplies signatures only for
            // its internal edges; routing its boundary would erase them.
            expression = expression.replace_multiple(graph.normal_emr_replacement(
                &region.remaining,
                &basis,
                &[W_.x___],
                |_| true,
            ));
        }
        let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        for (pair, edge, data) in graph.iter_edges_of(&region.remaining) {
            if !pair.is_paired() {
                continue;
            }
            ensure!(
                data.data.mass_atom(model).is_zero(),
                "R* reference preparation currently requires massless physical propagators"
            );
            let hard = basis.loop_atom::<Atom>(edge, GS.emr_mom, &[], true);
            ensure!(
                !hard.is_zero(),
                "an R* remaining propagator has no hard loop momentum"
            );
            let momentum =
                basis.loop_atom(edge, GS.emr_mom, std::slice::from_ref(&minkowski), true)
                    + basis.ext_atom(edge, GS.emr_mom, std::slice::from_ref(&minkowski), true);
            let mass_squared = if edge == marked {
                reference_mass_squared.clone()
            } else {
                Atom::Zero
            };
            // `den` differentiates only its fourth argument. The third is
            // converter metadata and is the mass itself, squared by Vakint's
            // existing graph bridge after the Taylor operation.
            expression /= function!(
                GS.den,
                usize::from(edge),
                hard,
                mass_squared.sqrt(),
                function!(SPENSO_TAG.dot, &momentum, &momentum) - &mass_squared
            );
        }
        let mut converter_graph = graph.clone();
        for slot in &region.insertions {
            let nodes = converter_graph
                .iter_nodes_of(&slot.support)
                .map(|(node, _, _)| node)
                .collect::<Vec<_>>();
            let first = nodes
                .first()
                .ok_or_else(|| eyre!("an R* contracted support contains no vertices"))?;
            let vertex = converter_graph.underlying[*first].clone();
            converter_graph.underlying.identify_nodes(&nodes, vertex);
        }
        converter_graph.underlying.forget_identification_history();
        drop(routing_stage);
        if let Some(input) = experiment_input {
            let routing = graph
                .normal_emr_replacement(&region.remaining, &basis, &[W_.x___], |_| true)
                .into_iter()
                .collect::<Vec<_>>();
            super::edge_experiment::compare(
                graph,
                &converter_graph,
                &region.remaining,
                &input,
                &routing,
                region.loops,
                &|view| Self::simplify_view(view),
            );
        }
        Ok(Self {
            region: region.clone(),
            routed: expression,
            contracted: OnceLock::new(),
            graph: converter_graph,
            basis,
            soft_edges,
        })
    }

    pub fn soft_edges(&self) -> &[EdgeIndex] {
        &self.soft_edges
    }

    /// The routed expression contracted to dots, as the outer and physical
    /// paths consume it. It is contracted once, on first use.
    pub fn expression(&self) -> Result<&Atom> {
        if let Some(contracted) = self.contracted.get() {
            return Ok(contracted);
        }
        let contracted = {
            let _algebra = StageSpan::new("gamma.remainder.algebra");
            Self::simplify(&self.routed)
                .wrap_err("while preparing the routed R* tensor expression")?
        };
        Ok(self.contracted.get_or_init(|| contracted))
    }

    /// Mutable access to the contracted expression, for tests decorating it.
    #[cfg(test)]
    pub fn expression_mut(&mut self) -> Result<&mut Atom> {
        self.expression()?;
        Ok(self.contracted.get_mut().expect("contracted above"))
    }

    /// Express a two-point amplitude in its first boundary momentum, keeping
    /// the edge orientation convention used by the graph's external routing.
    pub fn physical_external_routing(
        graph: &Graph,
        boundary: &SuBitGraph,
    ) -> Result<(EdgeIndex, Vec<Replacement>)> {
        let externals = graph
            .iter_edges_of(boundary)
            .map(|(pair, edge, _)| match pair {
                HedgePair::Unpaired { flow, .. } => Ok((edge, flow)),
                _ => Err(eyre!(
                    "physical pole reconstruction needs true external edges"
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        ensure!(
            externals.len() == 2,
            "physical pole reconstruction currently requires a two-point graph"
        );
        let (external, flow) = externals[0];
        let routing = externals
            .iter()
            .map(|(edge, other_flow)| {
                let sign = if *edge == external || *other_flow != flow {
                    1
                } else {
                    -1
                };
                Replacement::new(
                    function!(GS.emr_mom, usize::from(*edge), W_.x___).to_pattern(),
                    Atom::num(sign) * function!(GS.emr_mom, usize::from(external), W_.x___),
                )
            })
            .collect::<Vec<_>>();
        Ok((external, routing))
    }

    /// Integrate a prepared propagator graph at its original external momentum.
    ///
    /// The vacuum converter intentionally forgets external flows and merges
    /// equal-mass two-edge bonds. Here the explicit quotient routing supplies
    /// the complete propagator momenta directly, preserving the physical scale.
    /// Tensor reduction averages transverse to that external momentum. All
    /// insertion coefficients enter Vakint before its Laurent expansion, so a
    /// child pole automatically promotes the required scalar-integral depth.
    pub fn integrate_physical(
        &self,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        vakint_settings: &VakintSettings,
        unique_id: usize,
    ) -> Result<Atom> {
        self.integrate_physical_numerator(
            self.physical_numerator()?,
            vakint,
            uv_settings,
            vakint_settings,
            unique_id,
        )
    }

    pub(super) fn physical_numerator(&self) -> Result<Atom> {
        Ok(self
            .expression()?
            .replace(function!(GS.den, W_.x___))
            .with(Atom::one()))
    }

    /// Certify a closed numerator in one canonical denominator frame.
    ///
    /// Only the compact particle graph enters graph isomorphism. Its chosen
    /// map transports the complete numerator, which is compared as an exact
    /// Atom after typed dummy normalization. Remaining automorphisms can miss
    /// matches conservatively. Open or unsupported tensors form singletons;
    /// they never fall back to canonizing a large expression-incidence graph.
    pub(super) fn physical_identity(
        &self,
        model: &Model,
        numerator: &Atom,
    ) -> Result<Option<(crate::uv::counterterm_identity::CountertermKey, Atom)>> {
        if UNSUPPORTED_PHYSICAL_MOMENTA
            .iter()
            .any(|&momentum| numerator.contains_symbol(momentum))
        {
            return Ok(None);
        }
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let tensor = SymbolicTensor::infer(cooking.try_cook(numerator.as_view())?)?;
        if !tensor.structure().logical_slots().is_empty() {
            return Ok(None);
        }
        let mut shell = self.graph.clone();
        let edges = shell
            .iter_edges()
            .map(|(pair, edge, _)| (pair.is_paired(), edge))
            .collect::<Vec<_>>();
        for (internal, edge) in edges {
            let retained = self.region.contains_edge(&shell, edge);
            let data = &mut shell.underlying[edge];
            data.num = Autogen::explicit(Atom::one());
            data.dod = Autogen::explicit(0);
            if internal && !retained {
                data.is_dummy = true;
            }
        }
        let vertices = shell
            .underlying
            .iter_nodes()
            .map(|(vertex, _, _)| vertex)
            .collect::<Vec<_>>();
        for vertex in vertices {
            let data = &mut shell.underlying[vertex];
            data.num = Autogen::explicit(Atom::one());
            data.dod = Autogen::explicit(0);
            data.vertex_rule = None;
        }
        shell.global_prefactor.num = Atom::one();
        shell.global_prefactor.projector = Atom::one();
        shell.loop_momentum_basis = self.basis.clone();
        let identity =
            crate::uv::counterterm_identity::CountertermIdentity::for_fixed_boundary_skeleton(
                &shell,
                &self.region.remaining,
                model,
            )?;
        let mut canonical = identity.canonicalize_integrand(numerator)?.expand();
        let mut matcher = SlotMatcher::default();
        let mut indexed = false;
        canonical.visitor(&mut |part| {
            if matches!(matcher.classify(part), SlotMatch::Explicit(_)) {
                indexed = true;
                false
            } else {
                !indexed
            }
        });
        if indexed {
            let tensor = SymbolicTensor::infer(cooking.try_cook(canonical.as_view())?)?;
            // The shared typed canonicalizer retains sum scopes. Symbolica
            // canonizes each top-level summand separately and treats subtrees
            // without dummy indices as scalar leaves. Invalid tensor/interface
            // inputs remain errors rather than silently becoming singletons.
            let normalized = tensor.canonize()?;
            canonical = cooking.uncook(normalized.expression().as_view());
        }
        Ok(Some((identity.key, canonical)))
    }

    pub(super) fn integrate_physical_numerator(
        &self,
        mut numerator: Atom,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        vakint_settings: &VakintSettings,
        unique_id: usize,
    ) -> Result<Atom> {
        let _stage = StageSpan::new("gamma.physical_integral");
        let routing_stage = StageSpan::new("gamma.physical_integral.routing");
        let mut timing = IntegrationTimings::new();
        if self.region.loops == 0 {
            drop(routing_stage);
            timing.routing = timing.lap();
            let result = Self::finish_tensor(numerator, unique_id)?;
            timing.finish = timing.lap();
            timing.log(&self.graph, unique_id, 0, "physical");
            return Ok(result);
        }
        let (external, external_routing) =
            Self::physical_external_routing(&self.graph, &self.region.boundary)?;
        numerator = numerator.replace_multiple(&external_routing);
        let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let mut reducer = TensorReducer::new(Atom::var(GS.dim)).with_external_vector(function!(
            GS.emr_mom,
            usize::from(external),
            &minkowski
        ))?;
        for edge in &self.basis.loop_edges {
            reducer = reducer.with_integrated_vector(function!(
                GS.emr_mom,
                usize::from(*edge),
                &minkowski
            ))?;
        }
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        timing.routing = timing.lap();
        drop(routing_stage);
        let materialization_stage =
            StageSpan::new("gamma.physical_integral.tensor_materialization");
        // Compact slash momenta inside gamma chains still carry tensor rank.
        // Expose those vectors before projection, using the same materializer
        // as the vacuum bridge. Scalar products and their inverses stay compact.
        let explicit = Integrated::materialize_tensor(&numerator)?;
        let source = reducer.admit(cooking.try_cook(explicit.as_view())?)?;
        timing.materialization = timing.lap();
        drop(materialization_stage);
        let projection_stage = StageSpan::new("gamma.physical_integral.tensor_reduction");
        let reduced = reducer
            .reduce(&source)
            .wrap_err("could not reduce a physical propagator tensor")?
            .into_expression()?;
        timing.tensor_projection = timing.lap();
        drop(projection_stage);
        let algebra_stage = StageSpan::new("gamma.physical_integral.algebra");
        numerator = Self::finish_tensor(cooking.uncook(reduced.expression().as_view()), unique_id)?;
        // Tensor restoration can spell scalar products as g(p,k). Normalize
        // them before mapping to Vakint's dot notation; represented metrics
        // and inverse external invariants retain their tensor/scalar roles.
        let scalar_products =
            SymbolicTensor::infer_indexed(cooking.try_cook(numerator.as_view())?)?.to_dots()?;
        numerator = cooking.uncook(scalar_products.expression().as_view());
        timing.algebra = timing.lap();
        drop(algebra_stage);
        let conversion_stage = StageSpan::new("gamma.physical_integral.vakint_conversion");

        let mut routing = self
            .basis
            .loop_edges
            .iter()
            .enumerate()
            .map(|(index, edge)| {
                Replacement::new(
                    function!(GS.emr_mom, usize::from(*edge), W_.x___).to_pattern(),
                    function!(S.k, index + 1, W_.x___),
                )
            })
            .collect::<Vec<_>>();
        routing.push(Replacement::new(
            function!(GS.emr_mom, usize::from(external), W_.x___).to_pattern(),
            function!(S.p, usize::from(external), W_.x___),
        ));
        numerator = numerator
            .replace_multiple(&routing)
            .replace(function!(
                SPENSO_TAG.dot,
                function!(W_.a_, W_.a___, &minkowski),
                function!(W_.b_, W_.b___, &minkowski)
            ))
            .with(function!(
                S.dot,
                function!(W_.a_, W_.a___),
                function!(W_.b_, W_.b___)
            ))
            .replace(function!(
                ETS.metric,
                Minkowski {}.to_symbolic([W_.a__]),
                Minkowski {}.to_symbolic([W_.b__])
            ))
            .with(function!(
                S.metric,
                Minkowski {}.to_symbolic([W_.a__]),
                Minkowski {}.to_symbolic([W_.b__])
            ));
        let mut propagators = Atom::one();
        for (index, (pair, edge, _)) in self.graph.iter_edges_of(&self.region.remaining).enumerate()
        {
            let HedgePair::Paired { source, sink } = pair else {
                continue;
            };
            let momentum = (self.basis.loop_atom::<Atom>(edge, GS.emr_mom, &[], true)
                + self.basis.ext_atom::<Atom>(edge, GS.emr_mom, &[], true))
            .replace_multiple(&external_routing)
            .replace_multiple(&routing);
            propagators *= function!(
                S.prop,
                index + 1,
                function!(
                    S.edge,
                    usize::from(self.graph.node_id(source)),
                    usize::from(self.graph.node_id(sink))
                ),
                momentum,
                0,
                self.graph[edge].extra_data.vakint_edge_power.unwrap_or(1)
            );
        }
        numerator *=
            symbolica::parse!(&uv_settings.vakint.additional_normalization).pow(self.region.loops);
        let input = numerator * function!(S.topo, propagators);
        debug_tags!(#uv, #integrated, #rstar, #trace;
            stage = "rstar_physical_remainder",
            log.input = input,
            "Integrating a UV preparation at its physical external momentum"
        );
        let mut settings = vakint_settings.clone();
        settings.use_dot_product_notation = true;
        timing.conversion = timing.lap();
        drop(conversion_stage);
        let backend_stage = StageSpan::new("gamma.physical_integral.backend");
        let integrated = vakint
            .evaluate_integral(&settings, input.as_view())
            .wrap_err("while integrating a physical UV spinney remainder")?;
        timing.backend = timing.lap();
        drop(backend_stage);
        let result = Self::finish_tensor(integrated, unique_id)?;
        timing.finish = timing.lap();
        timing.log(&self.graph, unique_id, self.region.loops, "physical");
        Ok(result)
    }

    /// Put the bare hard term and its local counterterms in the same covariant
    /// basis before a quotient's pole projection. Hard scalar integrals remain
    /// unevaluated, including the dimension-dependent projector coefficients.
    pub fn project_hard(&self, expression: &Atom) -> Result<Atom> {
        project_hard_vacuum(expression, &self.basis.loop_edges.raw)
    }

    /// Hard Taylor expansion leaves the reference mass fixed. These momenta
    /// are external to this component, even when they are loops of its parent.
    ///
    /// The scaling and the truncation act on the routed product before any
    /// index is contracted, so only the kept Taylor terms are contracted
    /// (as `approx::local_4d` does for local counterterms). Contraction is
    /// linear and the parameter is a scalar, so this is the Taylor
    /// polynomial of the contracted expression, [`Self::taylor_contracted`].
    pub fn taylor(&self, degree: usize) -> Result<Atom> {
        let _stage = StageSpan::new("gamma.remainder.hard_taylor");
        let parameter = symbol!("rstar_hard_taylor"; Scalar);
        streamed::dump("hard_routed", self.routed.as_view());
        let scaled = Self::scale_edges(&self.routed, &self.soft_edges, parameter);
        let expanded = TaylorExpansion::new(parameter, degree).evaluate(&scaled)?;
        drop(scaled);
        let truncated = {
            let _algebra = StageSpan::new("gamma.remainder.hard_taylor_algebra");
            streamed::report_size("hard_taylor", expanded.as_view());
            streamed::dump("hard_taylor", expanded.as_view());
            Self::simplify_hard_taylor(&expanded)
                .wrap_err("while simplifying the R* hard Taylor polynomial")?
        };
        if cfg!(test) || hard_taylor_check() {
            let reference = self.taylor_contracted(degree)?;
            let identical = truncated == reference;
            let equal = identical || (&truncated - &reference).expand().is_zero();
            if cfg!(test) {
                assert!(
                    equal,
                    "hard Taylor before contraction differs from the contracted path\nnew: {truncated}\nold: {reference}"
                );
            } else {
                eprintln!(
                    "HARD_TAYLOR_CHECK loops={} degree={degree} identical={identical} equal={equal}",
                    self.region.loops
                );
            }
        }
        Ok(truncated)
    }

    /// [`Self::simplify`] of a hard Taylor polynomial. The factorized
    /// contraction of a deeply nested product of many index-sharing sums
    /// (27105: a 41 KB polynomial) can grow without bound; under an active
    /// resource budget (`BETA_PROCESS_RSS_LIMIT_GIB`) it aborts, and the
    /// polynomial is then expanded and simplified in chunks of summands:
    /// the same value, not the same factored form, so this path runs only
    /// where the whole simplification cannot finish
    /// (`GAMMALOOP_ALGEBRA_LEGACY=hard_taylor_fallback` lets the abort fail
    /// the diagram instead).
    fn simplify_hard_taylor(expanded: &Atom) -> Result<Atom> {
        let attempt =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| Self::simplify(expanded)));
        match attempt {
            Ok(result) => result,
            Err(payload) => {
                let message = payload
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| {
                        payload
                            .downcast_ref::<&str>()
                            .map(|text| (*text).to_owned())
                    })
                    .unwrap_or_default();
                if !message.starts_with("resource budget exhausted")
                    || streamed::legacy("hard_taylor_fallback")
                {
                    std::panic::resume_unwind(payload);
                }
                eprintln!(
                    "HARD_TAYLOR_FALLBACK bytes={} reason={message}",
                    expanded.as_view().get_byte_size()
                );
                // Return the aborted frontier's pages to the kernel before the
                // budget is checked again (the allocator keeps freed memory).
                if let Some(release) = MEMORY_RELEASE.get() {
                    release();
                    std::thread::sleep(std::time::Duration::from_millis(1500));
                }
                let started = std::time::Instant::now();
                let expanded = expanded.expand();
                let simplified = Self::simplify_summands(expanded.as_view())?;
                eprintln!(
                    "HARD_TAYLOR_FALLBACK done terms_bytes={} result_bytes={} seconds={:.1}",
                    expanded.as_view().get_byte_size(),
                    simplified.as_view().get_byte_size(),
                    started.elapsed().as_secs_f64()
                );
                Ok(simplified)
            }
        }
    }

    /// The previous hard Taylor path: contract the whole routed expression,
    /// then scale and truncate it. Kept as the reference for [`Self::taylor`].
    fn taylor_contracted(&self, degree: usize) -> Result<Atom> {
        let parameter = symbol!("rstar_hard_taylor"; Scalar);
        let scaled = Self::scale_edges(self.expression()?, &self.soft_edges, parameter);
        let expanded = TaylorExpansion::new(parameter, degree).evaluate(&scaled)?;
        Self::simplify(&expanded).wrap_err("while simplifying the R* hard Taylor polynomial")
    }

    #[allow(clippy::too_many_arguments)]
    pub fn integrate(
        &self,
        expression: &Atom,
        root_graph: &Graph,
        degree: usize,
        reference_mass_squared: Symbol,
        auxiliary_mass_squared: &Atom,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        vakint_settings: &VakintSettings,
        unique_id: usize,
        root_projector: Option<&RootProjector<'_>>,
    ) -> Result<Atom> {
        let mut timing = IntegrationTimings::new();
        let root_external = root_graph.external_filter::<SuBitGraph>();
        let root_edges = root_graph
            .iter_edges_of(&root_external)
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        // Quotient UV maps return original Q labels. Route newly introduced
        // component edges, preserving this component's independent hard loops.
        let routed = if self.region.loops == 0 {
            Cow::Borrowed(expression)
        } else {
            Cow::Owned(
                expression.replace_multiple(self.graph.normal_emr_replacement(
                    &self.region.remaining,
                    &self.basis,
                    &[W_.x___],
                    HedgePair::is_paired,
                )),
            )
        };
        let allowed_edges = self
            .basis
            .loop_edges
            .iter()
            .chain(&root_edges)
            .map(|edge| usize::from(*edge))
            .collect::<BTreeSet<_>>();
        let mut unresolved = BTreeSet::new();
        routed.visitor(&mut |part| {
            if let AtomView::Fun(momentum) = part
                && momentum.get_symbol() == GS.emr_mom
                && let Ok(edge) = usize::try_from(momentum.get(0))
                && !allowed_edges.contains(&edge)
            {
                unresolved.insert(edge);
            }
            true
        });
        ensure!(
            unresolved.is_empty(),
            "an R* quotient left unintegrated boundary loop carriers: {unresolved:?}"
        );
        timing.routing = timing.lap();
        let parameter = symbol!("rstar_outer_taylor"; Scalar);
        let taylor = TaylorExpansion::new(parameter, degree);
        // Each intermediate is dropped as soon as its successor exists: on the
        // largest diagrams each of them is hundreds of megabytes.
        let scaled = Self::scale_edges(&routed, &root_edges, parameter);
        drop(routed);
        let compensated =
            taylor.compensate_mass(&scaled, reference_mass_squared, auxiliary_mass_squared);
        drop(scaled);
        let polynomial = taylor.evaluate(&compensated)?;
        drop(compensated);
        timing.taylor = timing.lap();
        let expanded = Self::simplify_summands(polynomial.as_view())
            .wrap_err("while simplifying the compensated R* outer Taylor polynomial")?;
        drop(polynomial);
        timing.algebra = timing.lap();
        // The root is the only boundary which may be projected here. Child
        // and quotient pole maps are already complete, and no later pole map
        // can multiply finite contraction terms by another insertion pole.
        // Project after Taylor expansion so physical external momenta in a
        // projector are spectators, rather than differentiated source labels.
        let expanded = if let Some(projector) = root_projector {
            let projected = projector(&expanded)
                .wrap_err("while projecting the complete R* outer Taylor polynomial")?;
            drop(expanded);
            Integrated::dimensionally_regularized(&projected)
        } else {
            expanded
        };
        timing.tensor_projection = timing.lap();
        let integrated = if self.region.loops == 0 {
            expanded
        } else {
            // By value: to_vakint drops it once it has been materialized.
            let input = to_vakint_integrand(
                expanded,
                &self.graph.underlying,
                &self.region.remaining,
                &self.graph.empty_subgraph::<SuBitGraph>(),
                &uv_settings.vakint,
                false,
            )?;
            timing.conversion = timing.lap();
            let integrated = integrate_single_mass(vakint, vakint_settings, &input)?;
            timing.backend = timing.lap();
            integrated
        };
        let result = Self::finish_tensor(integrated, unique_id)?;
        timing.finish = timing.lap();
        timing.log(root_graph, unique_id, self.region.loops, "single_mass");
        Ok(result)
    }

    /// Use exactly the ordinary R operation's tensor and dimension restoration.
    pub fn finish_tensor(expression: Atom, unique_id: usize) -> Result<Atom> {
        let _stage = StageSpan::new("gamma.tensor_restoration");
        // Restore compact external scalar products before opening Vakint's
        // remaining dot shorthand. Both physical quotients and normalized
        // root projectors can contain reciprocals of these invariants;
        // expanding them would create inverse indexed tensors. Mixed loop
        // products and already indexed arguments retain the strict converter.
        let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let expression = expression
            .replace(function!(
                S.dot,
                function!(S.p, W_.i_),
                function!(S.p, W_.j_)
            ))
            .with(function!(
                SPENSO_TAG.dot,
                function!(S.p, W_.i_, &minkowski),
                function!(S.p, W_.j_, &minkowski)
            ));
        let indexed = Vakint::convert_from_dot_notation(expression.as_view());
        let restored = Integrated::restore_numerator(indexed, unique_id);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let reduced = SymbolicTensor::infer_indexed(cooking.try_cook(restored.as_view())?)?
            .contract(ContractSettings {
                rank_one: false,
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })?
            .to_dots()?;
        let regulated =
            Integrated::dimensionally_regularized(&cooking.uncook(reduced.expression().as_view()));
        // Keep collected trace/chain powers for their domain algebra. Opening
        // a powered trace first would give every copy the same spinor dummies,
        // turning two independent closed traces into repeated gamma factors.
        let collected = SymbolicTensor::infer_indexed(cooking.try_cook(regulated.as_view())?)?
            .contract(ContractSettings::default())?
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::canonical()),
                epsilon: true,
                color: Some(ColorSimplifySettings::default()),
                ..Default::default()
            })?;
        // Clifford contraction can create new scalar dimensions. Expand
        // those before the caller takes poles or inserts this tensor in an
        // enclosing graph; open representation dimensions remain metadata.
        let reduced = cooking.uncook(collected.expression().as_view());
        Integrated::warn_repeated_lorentz_indices(
            &reduced,
            "R*",
            Some(collected.reduction_status()),
        );
        Ok(Integrated::dimensionally_regularized(&reduced))
    }

    fn scale_edges(expression: &Atom, edges: &[EdgeIndex], parameter: Symbol) -> Atom {
        let replacements = edges
            .iter()
            .map(|edge| {
                let momentum = function!(GS.emr_mom, usize::from(*edge), W_.x___);
                Replacement::new(momentum.to_pattern(), Atom::var(parameter) * momentum)
            })
            .collect::<Vec<_>>();
        expression.replace_multiple(&replacements)
    }

    pub(super) fn simplify(expression: &Atom) -> Result<Atom> {
        debug_tags!(#uv, #integrated, #rstar, #trace;
            stage = "rstar_tensor_input",
            log.expression = expression,
            "Preparing an explicit R* expression for native tensor algebra"
        );
        Self::simplify_view(expression.as_view())
    }

    /// [`Self::simplify`] summand by summand: the simplification acts on
    /// every summand alone, and a normalized sum is canonical. A large sum
    /// is simplified in chunks of summands and added up once, unless a chunk
    /// introduces a dummy index, which makes the whole sum be simplified.
    fn simplify_summands(expression: AtomView<'_>) -> Result<Atom> {
        let limit = streamed::piece_limit();
        if expression.get_byte_size() < limit || !matches!(expression, AtomView::Add(_)) {
            return Self::simplify_view(expression);
        }
        let mut accumulated = streamed::SumAccumulator::default();
        let mut introduced = false;
        let outcome = streamed::for_each_summand_chunk(expression, limit, &mut |chunk| {
            let simplified = Self::simplify_view(chunk)?;
            if streamed::introduces_dummies(chunk, simplified.as_view()) {
                introduced = true;
                return Err(eyre!("a simplification introduced a dummy index"));
            }
            accumulated.push(simplified);
            Ok(())
        });
        let simplified = match outcome {
            Ok(()) => accumulated.finish(),
            Err(_) if introduced => return Self::simplify_view(expression),
            Err(error) => return Err(error),
        };
        if streamed::check_streamed(expression) {
            let whole = Self::simplify_view(expression)?;
            streamed::report_check("rstar_simplify", expression, &simplified, &whole);
        }
        Ok(simplified)
    }

    fn simplify_view(expression: AtomView<'_>) -> Result<Atom> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        // As in the ordinary R simplifier, retain open color tensors here.
        // Reduce color identities after the quotient has closed their indices.
        // Leave trace/chain powers intact for the same reason as finish_tensor:
        // their internal spinor sums must remain independent between copies.
        let simplified = SymbolicTensor::infer_indexed(cooking.try_cook(expression)?)?
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::canonical()),
                epsilon: true,
                ..Default::default()
            })?
            .to_dots()?;
        Ok(cooking.uncook(simplified.expression().as_view()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        finalized_runtime_dot, graph::parse::IntoFinalizedRuntimeGraph,
        initialisation::test_initialise, utils::load_generic_model, uv::rstar::RStarPlan,
    };
    use vakint::symbols::S;

    #[test]
    fn simplify_powered_dirac_trace_preserves_open_color_tensor() -> Result<()> {
        use crate::numerator::aind::Aind;
        use idenso::{
            IndexTooling,
            representations::{Bispinor, ColorAdjoint, ColorAntiFundamental, ColorFundamental},
        };
        use spenso::structure::slot::{DummyAind, ParseableAind};

        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let trace = spenso::trace!(
            Bispinor {}.new_rep(4).to_symbolic([]),
            idenso::gamma!(&mu),
            idenso::gamma!(&nu)
        );
        let spectator = idenso::color_t!(
            ColorAdjoint {}
                .new_rep(8)
                .to_symbolic([Aind::Hedge(70, 1).to_atom()]),
            ColorFundamental {}
                .new_rep(3)
                .to_symbolic([Aind::Hedge(71, 1).to_atom()]),
            ColorAntiFundamental {}
                .new_rep(3)
                .to_symbolic([Aind::Hedge(72, 1).to_atom()])
        );
        let source = trace.pow(2) * &spectator;
        let reduced = RStarIntegrand::simplify(&source)?;
        assert!(!reduced.has_repeated_lorentz_indices());
        // The early Taylor boundary retains symbolic D and the open color
        // generator; only its two independent Dirac traces are evaluated.
        let expected = Atom::num(16) * Atom::var(GS.dim) * spectator;
        super::super::tests::assert_same_counterterm(&reduced, &expected)?;
        Ok(())
    }

    #[test]
    fn finish_tensor_two_traces_contract_shared_lorentz_indices_before_poles() -> Result<()> {
        use crate::numerator::aind::Aind;
        use idenso::{IndexTooling, representations::Bispinor};
        use spenso::structure::slot::{DummyAind, ParseableAind};

        test_initialise()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let ports: [_; 4] =
            std::array::from_fn(|i| spin.to_symbolic([Aind::Hedge((i + 3) as u16, 1).to_atom()]));
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let gamma =
            |a: &Atom, b: &Atom, index: &Atom| function!(idenso::dirac::AGS.gamma, a, b, index);
        let source = gamma(&ports[0], &ports[1], &mu)
            * gamma(&ports[1], &ports[0], &nu)
            * gamma(&ports[2], &ports[3], &mu)
            * gamma(&ports[3], &ports[2], &nu)
            / epsilon.pow(2);
        let reduced = RStarIntegrand::finish_tensor(source, 9954)?;
        assert!(!reduced.has_repeated_lorentz_indices());
        let poles = reduced.series(GS.dim_epsilon, 0, -1)?.to_atom();
        let expected = Atom::num(64) / epsilon.pow(2) - Atom::num(32) / epsilon;
        super::super::tests::assert_same_counterterm(&poles, &expected)?;
        Ok(())
    }

    #[test]
    fn finish_tensor_reduces_clifford_words_before_caller_poles() -> Result<()> {
        use crate::numerator::aind::Aind;
        use idenso::representations::Bispinor;
        use spenso::structure::slot::{DummyAind, ParseableAind};

        test_initialise()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let a = spin.to_symbolic([Aind::Hedge(3, 1).to_atom()]);
        let b = spin.to_symbolic([Aind::Hedge(4, 1).to_atom()]);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let momentum = function!(GS.emr_mom, 7, mink.to_symbolic([]));
        let gamma = |argument: &Atom| {
            function!(
                idenso::dirac::AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                argument
            )
        };
        let word = spenso::chain!(&a, &b, gamma(&mu), gamma(&momentum), gamma(&mu));
        let source = word / epsilon.pow(2);
        let actual = RStarIntegrand::finish_tensor(source, 9951)?
            .series(GS.dim_epsilon, 0, -1)?
            .to_atom();
        let expected = (Atom::num(-2) + Atom::num(2) * &epsilon)
            * function!(idenso::dirac::AGS.gamma, &a, &b, &momentum)
            / epsilon.pow(2);
        super::super::tests::assert_same_counterterm(&actual, &expected)?;
        Ok(())
    }

    #[test]
    fn finish_tensor_preserves_external_scalar_denominators() -> Result<()> {
        use crate::uv::rstar::tests::assert_same_counterterm;

        test_initialise()?;
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let mu = minkowski.to_symbolic([Atom::var(symbol!("rstar_scalar_exit::mu"))]);
        let nu = minkowski.to_symbolic([Atom::var(symbol!("rstar_scalar_exit::nu"))]);
        let spectator = function!(ETS.metric, &mu, &nu);
        let source = |a, b| function!(S.dot, function!(S.p, a), function!(S.p, b));
        let restored = |a, b| {
            function!(
                SPENSO_TAG.dot,
                function!(GS.emr_mom, a, minkowski.to_symbolic([])),
                function!(GS.emr_mom, b, minkowski.to_symbolic([]))
            )
        };
        // Cover a diagonal denominator, a mixed external invariant and a
        // reciprocal polynomial, with a retained open tensor boundary.
        for (input, expected) in [
            (source(0, 0).pow(-1), restored(0, 0).pow(-1)),
            (source(0, 1).pow(-2), restored(0, 1).pow(-2)),
            (
                (source(0, 0) + source(0, 1)).pow(-2),
                (restored(0, 0) + restored(0, 1)).pow(-2),
            ),
        ] {
            let actual = RStarIntegrand::finish_tensor(input * &spectator, 9911)?;
            assert_same_counterterm(&actual, &(expected * &spectator))?;
        }
        // Scalar restoration must not conceal unsupported inverse tensors
        // or an inverse product carrying unintegrated loop momentum.
        for invalid in [
            function!(S.p, 0, 77).pow(-1),
            function!(S.dot, function!(S.k, 1), function!(S.p, 0)).pow(-1),
        ] {
            assert!(RStarIntegrand::finish_tensor(invalid, 9912).is_err());
        }
        Ok(())
    }

    #[test]
    fn quadratic_gluon_hard_preparation_has_no_uv_pole_at_fixed_mass() -> Result<()> {
        use crate::uv::rstar::tests::{gluon_self_energy, parallel_pair, settings};

        test_initialise()?;
        let model = load_generic_model("sm");
        let graph =
            gluon_self_energy(&model, 2, 4, "V_36", |graph| parallel_pair(graph).is_some())?;
        let bubble = parallel_pair(&graph).unwrap();
        let settings = settings();
        let vakint = Vakint::new()?;
        let marked = bubble[0];
        let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
        let correction = plan
            .ir
            .iter()
            .find(|correction| correction.taylor_order == 2 && graph.n_loops(&correction.hard) == 1)
            .ok_or_else(|| eyre!("fixture has no quadratic one-loop hard component"))?;
        let bare = correction
            .preparation
            .iter()
            .find(|term| term.counterterms.is_empty())
            .ok_or_else(|| eyre!("hard preparation is missing its bare term"))?;
        let mass_squared = symbol!("rstar_fixed_hard_mass_squared"; Scalar);
        let children = plan.evaluate_uv_children(
            &graph,
            &model,
            &Atom::var(mass_squared).sqrt(),
            &vakint,
            &settings,
        )?;
        let hard = RStarIntegrand::new(
            &graph,
            &model,
            &bare.remainder,
            &[],
            marked,
            Atom::var(mass_squared),
        )?;
        let polynomial = hard.taylor(correction.taylor_order as usize)?;
        let input = to_vakint_integrand(
            &polynomial,
            &hard.graph.underlying,
            &hard.region.remaining,
            &hard.graph.empty_subgraph::<SuBitGraph>(),
            &settings.vakint,
            false,
        )?;
        let mut integration_settings = settings.vakint.true_settings();
        integration_settings.number_of_terms_in_epsilon_expansion = 4;
        let integrated = RStarIntegrand::finish_tensor(
            integrate_single_mass(&vakint, &integration_settings, &input)?,
            9851,
        )?;
        let counterterm = children
            .get(&correction.hard)
            .ok_or_else(|| eyre!("hard component has no completed massive UV counterterm"))?;
        let prepared = RStarIntegrand::finish_tensor(&integrated + counterterm, 9852)?;
        let poles = prepared.series(GS.dim_epsilon, 0, -1)?.to_atom();
        debug_tags!(#uv, #integrated, #rstar, #trace;
            stage = "rstar_fixed_mass_hard_preparation",
            log.bare = integrated,
            log.counterterm = counterterm,
            log.poles = poles,
            "Checking hard preparation before any quotient UV map or outer Taylor operation"
        );
        assert!(
            poles.together().is_zero(),
            "fixed-mass hard preparation retained UV poles: {poles}"
        );
        Ok(())
    }

    /// Hard Taylor before contraction against the contracted reference, on
    /// every hard preparation term of gluonic two- and three-loop graphs,
    /// bare and with completed counterterm insertions, at every order up to
    /// the one the plan uses.
    #[test]
    fn hard_taylor_before_contraction_matches_the_contracted_path() -> Result<()> {
        use crate::uv::rstar::tests::{gluon_self_energy, parallel_pair, settings};

        test_initialise()?;
        let model = load_generic_model("sm");
        let settings = settings();
        let vakint = Vakint::new()?;
        let mass_squared = symbol!("rstar_hard_check_mass_squared"; Scalar);
        let graphs = [
            gluon_self_energy(&model, 2, 4, "V_36", |graph| parallel_pair(graph).is_some())?,
            gluon_self_energy(&model, 3, 6, "V_36", |graph| parallel_pair(graph).is_some())?,
        ];
        let mut compared = 0;
        let mut with_insertions = 0;
        for graph in &graphs {
            let marked = parallel_pair(graph).unwrap()[0];
            let plan = RStarPlan::new(graph, &model, marked, &settings)?;
            let children = plan.evaluate_uv_children(
                graph,
                &model,
                &Atom::var(mass_squared).sqrt(),
                &vakint,
                &settings,
            )?;
            for correction in &plan.ir {
                for term in &correction.preparation {
                    if term.remainder.is_scaleless(graph, plan.marked_edge) {
                        continue;
                    }
                    let Some(insertions) = term
                        .counterterms
                        .iter()
                        .map(|child| children.get(&child.argument.support).cloned())
                        .collect::<Option<Vec<_>>>()
                    else {
                        continue;
                    };
                    let hard = RStarIntegrand::new(
                        graph,
                        &model,
                        &term.remainder,
                        &insertions,
                        plan.marked_edge,
                        Atom::var(mass_squared),
                    )?;
                    for degree in 0..=correction.taylor_order as usize {
                        let new = hard.taylor(degree)?;
                        let old = hard.taylor_contracted(degree)?;
                        assert!(
                            (&new - &old).expand().is_zero(),
                            "hard Taylor order {degree} differs:\nnew: {new}\nold: {old}"
                        );
                        compared += 1;
                    }
                    with_insertions += usize::from(!insertions.is_empty());
                }
            }
        }
        assert!(
            compared > 0 && with_insertions > 0,
            "{compared} {with_insertions}"
        );
        Ok(())
    }

    #[test]
    fn hard_taylor_preserves_the_single_mass_squared_in_the_existing_bridge() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph: Graph = finalized_runtime_dot!(
            digraph bubble {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=2 sink="{ufo_order:0}"];
                B:1 -> e [id=3 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, &model
        )?;
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            graph.underlying[edge].particle =
                graph[edge].particle.clone().override_mass(Some(Atom::Zero));
        }
        let settings = UVgenerationSettings::default();
        let marked = EdgeIndex::from(0);
        let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
        let mass_squared = symbol!("rstar_bridge_mass_squared");
        let input = RStarIntegrand::new(
            &graph,
            &model,
            &plan.ordinary_uv[0].remainder,
            &[],
            marked,
            Atom::var(mass_squared),
        )?;
        assert!(!input.soft_edges().is_empty());
        // A positive Taylor order differentiates shifted Lorentz products.
        // Its parameter must be scalar so it factors out of those products
        // before Symbolica computes the series.
        let quadratic = input.taylor(2)?;
        assert!(input.soft_edges().iter().any(|edge| {
            quadratic
                .replace(function!(GS.emr_mom, usize::from(*edge), W_.x___))
                .match_iter()
                .next()
                .is_some()
        }));
        let hard = input.taylor(0)?;
        for edge in input.soft_edges() {
            assert!(
                hard.replace(function!(GS.emr_mom, usize::from(*edge), W_.x___))
                    .match_iter()
                    .next()
                    .is_none()
            );
        }
        let converted = to_vakint_integrand(
            &hard,
            &input.graph.underlying,
            &input.region.remaining,
            &input.graph.empty_subgraph::<SuBitGraph>(),
            &settings.vakint,
            false,
        )?;
        assert_eq!(converted.0.len(), 1);
        let mut masses = BTreeSet::new();
        converted.0[0].integral.visitor(&mut |part| {
            if let AtomView::Fun(prop) = part
                && prop.get_symbol() == S.prop
            {
                masses.insert(prop.get(3).to_owned());
                return false;
            }
            true
        });
        assert_eq!(
            masses,
            BTreeSet::from([Atom::Zero, Atom::var(mass_squared)])
        );
        // Boundary sources are retained individually, including the edge that
        // is dependent in the internal momentum routing. They must never be
        // replaced by the zero signatures of an internal-subgraph basis.
        for edge in input.soft_edges() {
            let vector = function!(
                GS.emr_mom,
                usize::from(*edge),
                Minkowski {}.new_rep(GS.dim).to_symbolic([])
            );
            let factor = function!(SPENSO_TAG.dot, &vector, &vector);
            let mut decorated = graph.clone();
            decorated.underlying[marked].num = crate::graph::Autogen::explicit(factor.clone());
            let preserved = RStarIntegrand::new(
                &decorated,
                &model,
                &input.region,
                &[],
                marked,
                Atom::var(mass_squared),
            )?;
            assert!(
                (preserved.expression()? - input.expression()? * &factor)
                    .together()
                    .is_zero()
            );
        }
        Ok(())
    }

    #[test]
    fn tensor_finish_restores_dots_and_regulates_only_scalar_dimension() -> Result<()> {
        test_initialise()?;
        let raw = function!(S.dot, function!(S.p, 7), function!(S.p, 7)) * function!(S.g, 1, 1);
        let actual = RStarIntegrand::simplify(&RStarIntegrand::finish_tensor(raw, 543)?)?;
        let vector = function!(GS.emr_mom, 7, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
        let expected = (Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon))
            * function!(SPENSO_TAG.dot, &vector, &vector);
        assert!(
            (actual - expected).expand().is_zero(),
            "Vakint spectators must retain the original native Lorentz representation"
        );
        Ok(())
    }

    #[test]
    fn physical_remainder_preserves_scale_and_promotes_insertion_poles() -> Result<()> {
        use crate::numerator::aind::Aind;
        use spenso::structure::slot::ParseableAind;

        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph: Graph = finalized_runtime_dot!(
            digraph bubble {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=2 sink="{ufo_order:0}"];
                B:1 -> e [id=3 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, &model
        )?;
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            graph.underlying[edge].particle =
                graph[edge].particle.clone().override_mass(Some(Atom::Zero));
        }
        let settings = super::super::tests::settings();
        let marked = EdgeIndex::from(0);
        let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
        let mut input = RStarIntegrand::new(
            &graph,
            &model,
            &plan.ordinary_uv[0].remainder,
            &[],
            marked,
            Atom::Zero,
        )?;
        let epsilon = Atom::var(GS.dim_epsilon);
        *input.expression_mut()? /= epsilon.pow(2);
        let mut integration = settings.vakint.true_settings();
        integration.integral_normalization_factor = vakint::LoopNormalizationFactor::GScheme;
        integration.number_of_terms_in_epsilon_expansion = 1;
        let vakint = Vakint::new()?;
        let actual = input.integrate_physical(&vakint, &settings, &integration, 5431)?;
        let actual = super::super::tests::canonical_tensor(&actual)?;
        let external = function!(
            GS.emr_mom,
            usize::from(input.soft_edges[0]),
            Minkowski {}.new_rep(GS.dim).to_symbolic([])
        );
        let logarithm = function!(
            Symbol::LOG,
            -function!(SPENSO_TAG.dot, &external, &external) / Atom::var(GS.mu_r_sq)
        );
        let expected = epsilon.pow(-3) - &logarithm / epsilon.pow(2)
            + logarithm.pow(2) / (Atom::num(2) * &epsilon);
        assert!(
            (&actual - &expected).expand().is_zero(),
            "physical remainders must retain scale logarithms and positive orders promoted by poles: actual={actual}, expected={expected}"
        );
        // A scalar loop-external product can be restored as a compact metric
        // by finish_tensor. It must reach Forcer as a dot, while an open
        // Lorentz metric and a reciprocal external invariant stay spectators.
        let left = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::Hedge(87, 0).to_atom()]);
        let right = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::Hedge(88, 0).to_atom()]);
        let metric = function!(ETS.metric, &left, &right);
        let loop_momentum = function!(
            GS.emr_mom,
            usize::from(input.basis.loop_edges.raw[0]),
            Minkowski {}.new_rep(GS.dim).to_symbolic([])
        );
        let original = input.expression()?.clone();
        *input.expression_mut()? *= function!(SPENSO_TAG.dot, &loop_momentum, &external)
            / function!(SPENSO_TAG.dot, &external, &external)
            * &metric;
        let actual_scalar = input.integrate_physical(&vakint, &settings, &integration, 5433)?;
        let actual_scalar = super::super::tests::canonical_tensor(&actual_scalar)?;
        let expected_scalar = &expected * metric / Atom::num(2);
        assert!(
            (&actual_scalar - &expected_scalar).expand().is_zero(),
            "physical scalar products and their tensor spectators must survive ingress: actual={actual_scalar}, expected={expected_scalar}"
        );
        *input.expression_mut()? = original;
        // A vacuum projector would wrongly discard this odd loop tensor.
        // With the physical external momentum, symmetry of the two lines
        // instead gives integral k^mu = p^mu/2 times the scalar bubble.
        let slot = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::Hedge(86, 0).to_atom()]);
        let loop_vector = function!(
            GS.emr_mom,
            usize::from(input.basis.loop_edges.raw[0]),
            &slot
        );
        *input.expression_mut()? *= loop_vector;
        let actual_tensor = input.integrate_physical(&vakint, &settings, &integration, 5432)?;
        let actual_tensor = super::super::tests::canonical_tensor(&actual_tensor)?;
        let expected_tensor = expected
            * function!(GS.emr_mom, usize::from(input.soft_edges[0]), &slot)
            / Atom::num(2);
        assert!(
            (&actual_tensor - &expected_tensor).expand().is_zero(),
            "the physical tensor projector must retain the external direction: actual={actual_tensor}, expected={expected_tensor}"
        );
        Ok(())
    }
}
