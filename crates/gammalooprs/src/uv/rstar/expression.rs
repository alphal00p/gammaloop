//! Explicit hard Taylor coefficients with the existing UV tensor conventions.

use std::collections::BTreeSet;

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
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::Replacement,
    symbol,
};
use vakint::{Vakint, VakintSettings, symbols::S};

use crate::{
    debug_tags,
    graph::{Graph, LMBext, LoopMomentumBasis},
    model::Model,
    utils::{GS, W_},
    uv::{
        UVgenerationSettings, UltravioletGraph,
        approx::integrated::{Integrated, to_vakint_integrand},
        uv_graph::UVE,
    },
};

use super::{GraphRegion, TaylorExpansion, integrate_single_mass, project::project_hard_vacuum};

/// The original tensor slots and loop carriers survive all graph contractions.
/// Only propagator endpoints are identified in the converter's graph.
pub(super) struct RStarIntegrand {
    pub region: GraphRegion,
    pub expression: Atom,
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
        let mut expression = graph
            .numerator(&region.remaining, &contracted)
            .to_d_dim(GS.dim)
            .get_single_atom()?;
        for insertion in insertions {
            expression *= insertion;
        }
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
        Ok(Self {
            region: region.clone(),
            expression: Self::simplify(&expression)
                .wrap_err("while preparing the routed R* tensor expression")?,
            graph: converter_graph,
            basis,
            soft_edges,
        })
    }

    pub fn soft_edges(&self) -> &[EdgeIndex] {
        &self.soft_edges
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
        let mut numerator = self
            .expression
            .replace(function!(GS.den, W_.x___))
            .with(Atom::one());
        if self.region.loops == 0 {
            return Self::finish_tensor(numerator, unique_id);
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
        let source = reducer.admit(cooking.try_cook(numerator.as_view())?)?;
        let reduced = reducer
            .reduce(&source)
            .wrap_err("could not reduce a physical propagator tensor")?
            .into_expression()?;
        numerator = Self::finish_tensor(cooking.uncook(reduced.expression().as_view()), unique_id)?;

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
                1
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
        let integrated = vakint
            .evaluate_integral(&settings, input.as_view())
            .wrap_err("while integrating a physical UV spinney remainder")?;
        // Physical quotients can contain inverse powers of p². Preserve this
        // scalar invariant: expanding it into indexed vectors before taking
        // its reciprocal would look like a negative power of an open tensor.
        let external_head = function!(S.p, usize::from(external));
        let external_vector = function!(S.p, usize::from(external), &minkowski);
        let integrated = integrated
            .replace(function!(S.dot, &external_head, &external_head))
            .with(function!(
                SPENSO_TAG.dot,
                &external_vector,
                &external_vector
            ));
        Self::finish_tensor(integrated, unique_id)
    }

    /// Put the bare hard term and its local counterterms in the same covariant
    /// basis before a quotient's pole projection. Hard scalar integrals remain
    /// unevaluated, including the dimension-dependent projector coefficients.
    pub fn project_hard(&self, expression: &Atom) -> Result<Atom> {
        project_hard_vacuum(expression, &self.basis.loop_edges.raw)
    }

    /// Hard Taylor expansion leaves the reference mass fixed. These momenta
    /// are external to this component, even when they are loops of its parent.
    pub fn taylor(&self, degree: usize) -> Result<Atom> {
        let parameter = symbol!("rstar_hard_taylor"; Scalar);
        let scaled = Self::scale_edges(&self.expression, &self.soft_edges, parameter);
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
    ) -> Result<Atom> {
        let root_external = root_graph.external_filter::<SuBitGraph>();
        let root_edges = root_graph
            .iter_edges_of(&root_external)
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        // Quotient UV maps return original Q labels. Route newly introduced
        // component edges, preserving this component's independent hard loops.
        let routed = if self.region.loops == 0 {
            expression.clone()
        } else {
            expression.replace_multiple(self.graph.normal_emr_replacement(
                &self.region.remaining,
                &self.basis,
                &[W_.x___],
                HedgePair::is_paired,
            ))
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
        let parameter = symbol!("rstar_outer_taylor"; Scalar);
        let taylor = TaylorExpansion::new(parameter, degree);
        let scaled = Self::scale_edges(&routed, &root_edges, parameter);
        let compensated =
            taylor.compensate_mass(&scaled, reference_mass_squared, auxiliary_mass_squared);
        let expanded = Self::simplify(&taylor.evaluate(&compensated)?)
            .wrap_err("while simplifying the compensated R* outer Taylor polynomial")?;
        if self.region.loops == 0 {
            return Self::finish_tensor(expanded, unique_id);
        }
        let input = to_vakint_integrand(
            &expanded,
            &self.graph.underlying,
            &self.region.remaining,
            &self.graph.empty_subgraph::<SuBitGraph>(),
            &uv_settings.vakint,
            false,
        )?;
        let integrated = integrate_single_mass(vakint, vakint_settings, &input)?;
        Self::finish_tensor(integrated, unique_id)
    }

    /// Use exactly the ordinary R operation's tensor and dimension restoration.
    pub fn finish_tensor(expression: Atom, unique_id: usize) -> Result<Atom> {
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
        let collected = SymbolicTensor::infer_indexed(cooking.try_cook(regulated.as_view())?)?
            .contract(ContractSettings::default())?
            .undo_trace()?
            .undo_chain()?
            .simplify_algebra(&AlgebraSettings {
                color: Some(ColorSimplifySettings::default()),
                ..Default::default()
            })?;
        Ok(cooking.uncook(collected.expression().as_view()))
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

    fn simplify(expression: &Atom) -> Result<Atom> {
        debug_tags!(#uv, #integrated, #rstar, #trace;
            stage = "rstar_tensor_input",
            log.expression = expression,
            "Preparing an explicit R* expression for native tensor algebra"
        );
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        // As in the ordinary R simplifier, retain open color tensors here.
        // Reduce color identities after the quotient has closed their indices.
        let simplified = SymbolicTensor::infer_indexed(cooking.try_cook(expression.as_view())?)?
            .undo_trace()?
            .undo_chain()?
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
                (preserved.expression - &input.expression * factor)
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
        input.expression /= epsilon.pow(2);
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
        // A vacuum projector would wrongly discard this odd loop tensor.
        // With the physical external momentum, symmetry of the two lines
        // instead gives integral k^mu = p^mu/2 times the scalar bubble.
        let slot = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::Hedge(86, 0).to_atom()]);
        input.expression *= function!(
            GS.emr_mom,
            usize::from(input.basis.loop_edges.raw[0]),
            &slot
        );
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
