use feynkit_graph::DOD;
use std::ops::Deref;

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::{
    HedgeGraph, PowersetIterator,
    involution::EdgeIndex,
    subgraph::{InternalSubGraph, SuBitGraph, SubGraphLike, SubSetLike, SubSetOps, subset::SubSet},
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    domains::atom::AtomField,
    function,
    poly::series::Series,
};
use tracing::debug;

use crate::{
    graph::{Edge, FeynmanGraph, Graph, HedgeData, LMBext, LmbError, LoopMomentumBasis, Vertex},
    integrands::process::param_builder::ParamBuilderGraph,
    model::Model,
    momentum::sample::LoopIndex,
    numerator::{AppliedFeynmanRule, Numerator},
    settings::global::MediumSettings,
    utils::{GS, W_},
    uv::{UVgenerationSettings, UvGraph},
};

use super::{Spinney, Wood};

/// GammaLoop's UV operations on a graph, beyond the combinatorics of [`UvGraph`].
pub trait UltravioletGraph:
    UvGraph<Basis = LoopMomentumBasis> + LMBext + FeynmanGraph + ParamBuilderGraph
{
    ///Get the numerator of the graph.
    /// If multiply_prefactor is true, the numerator is multiplied by the global  prefactor. (just num not projector)
    fn numerator<S: SubGraphLike + SubSetOps>(
        &self,
        subgraph: &S,
        without: &S,
    ) -> Result<Numerator<AppliedFeynmanRule>>;
    fn denominator<S: SubGraphLike, T: Fn(&Edge) -> isize>(
        &self,
        subgraph: &S,
        model: &Model,
        edge_powers: T,
    ) -> Atom;

    /// The classified spinneys of `subgraph`, or only the empty spinney when UV
    /// subtraction is disabled.
    fn classified_spinneys<S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        model: &Model,
        settings: &UVgenerationSettings,
        lmb: &LoopMomentumBasis,
    ) -> Result<Vec<Spinney>> {
        if !settings.subtract_uv {
            return Ok(vec![Spinney::empty(self)]);
        }
        Ok(self.classify_spinneys(subgraph, model, &settings.renormalization_prescription, lmb)?)
    }

    fn remove_full_observable_spinney(&self, spinneys: &mut Vec<Spinney>, medium: &MediumSettings) {
        if !medium.vacuum_subtraction {
            return;
        }

        // Completed thermal weights already carry (1 - V). An operation covering
        // all loops leaves no medium dependence, even when tree edges survive.
        // Omit it before either 4D or 3D computation; spinneys contain no bridges.
        let graph = self.hedge_graph();
        let full_observable =
            InternalSubGraph::cleaned_filter_pessimist(graph.non_bridges(), graph);
        spinneys.retain(|spinney| {
            spinney.subgraph().is_empty() || *spinney.subgraph() != full_observable
        });
    }

    fn all_limits<E, V, H, S: SubGraphLike>(
        &self,
        model: &crate::model::Model,
        subgraph: &S,
        expr: &Atom,
        expansion: Symbol,
        lmb: &LoopMomentumBasis,
    ) -> Vec<(SubSet<LoopIndex>, Series<AtomField>)>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let mom_reps = self.normal_emr_replacement(subgraph, lmb, &[W_.x___], |_s| true);

        let ose_reps = self.get_ose_replacements(model);
        // for x in &mom_reps {
        //     println!("LMB replacement: {x}");
        // }

        // for r in &ose_reps {
        //     println!("ose{r}");
        // }

        // This separate analytic copy uses native square-root Laurent series:
        // after UV rescaling an energy invariant can have a pole, beyond the
        // compact function's ordinary Taylor derivative domain. Lower energies
        // after OSE substitution without distributing any numerator factors.
        let expr = expr
            .replace_multiple(&ose_reps)
            .replace_multiple(&mom_reps)
            .replace(function!(GS.broadcasting_sqrt, W_.a_))
            .with(Atom::var(W_.a_).sqrt())
            .replace(function!(GS.on_shell_energy, W_.a_, W_.prop_))
            .with(Atom::var(W_.prop_).sqrt());
        // .replace_multiple(&q3_reps);
        let mut loops = PowersetIterator::<LoopIndex>::new(lmb.loop_edges.len() as u8);

        let mut limits = Vec::new();

        loops.next();

        for ls in loops {
            let mut expr = expr.clone();
            for l in ls.included_iter() {
                let e = usize::from(lmb.loop_edges[l]) as i64;
                expr = expr
                    .replace(function!(GS.emr_mom, e, W_.x___))
                    .with(function!(GS.emr_mom, e, W_.x___) / expansion);

                expr /= Atom::var(expansion).pow(3);
            }

            let series = expr.series(expansion, Atom::Zero, 0).unwrap();

            // expr = series.to_atom().expand();

            limits.push((ls, series));
        }
        limits
    }

    fn wood<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        model: &Model,
    ) -> Result<Wood>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        self.wood_with_settings(
            subgraph,
            model,
            &UVgenerationSettings::default(),
            &self.as_ref().lmb_of(&self.as_ref().full_filter()),
        )
    }

    fn wood_with_settings<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        model: &Model,
        settings: &UVgenerationSettings,
        lmb: &LoopMomentumBasis,
    ) -> Result<Wood>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        Ok(Wood::from_spinneys(
            self.classified_spinneys(subgraph, model, settings, lmb)?,
            self,
        ))
    }

    fn compute_dod<S: SubGraphLike<Base = SuBitGraph> + SubSetOps>(&self, subgraph: &S) -> i32;

    /// Degree of divergence of the simultaneous loop-energy integrations.
    ///
    /// The existing UV DOD already contains the numerator and denominator
    /// scaling. Replacing the four-dimensional measure by one energy measure
    /// per loop therefore amounts to subtracting three powers per loop.
    fn compute_energy_dod<E, V, H, S>(&self, subgraph: &S) -> i32
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
        S: SubGraphLike<Base = SuBitGraph> + SubSetOps,
    {
        self.compute_dod(subgraph) - 3 * self.n_loops(subgraph) as i32
    }
}

impl AsRef<HedgeGraph<Edge, Vertex, HedgeData>> for Graph {
    fn as_ref(&self) -> &HedgeGraph<Edge, Vertex, HedgeData> {
        &self.underlying
    }
}

impl Graph {
    /// Reject source integrals whose loop-energy contours do not vanish at
    /// infinity. Check the existing unpruned cycle-union inventory so this
    /// source-level contract does not depend on UV-spinney pruning.
    pub(crate) fn ensure_energy_convergent_cycles<S>(&self, domain: &S) -> Result<()>
    where
        S: SubGraphLike<Base = SuBitGraph>,
    {
        let numerator = self.production_numerator_atom_for_full_3d_expression();
        if numerator.is_zero() {
            return Ok(());
        }
        let mut failures = self
            .all_cycle_unions(domain)
            .into_iter()
            .filter(|cycle| !cycle.is_empty())
            .filter_map(|cycle| {
                let loop_count = self.n_loops(&cycle.filter);
                // `compute_dod` deliberately covers only the edge and vertex
                // rules used by UV classification. Source convergence instead
                // sees the complete factorized production numerator, including
                // graph-global factors and momentum dependence stored outside
                // this cycle. The production helper excludes initial-cut
                // factors without changing UV DOD semantics or expanding the
                // numerator.
                let lmb = self.lmb_of(&cycle.filter);
                let integrand = &numerator / self.denominator(&cycle.filter, &self.model, |_| 1);
                let four_d_dod = self
                    .uv_rescaled(cycle.filter.included(), loop_count, &lmb, &lmb, &integrand)
                    .trailing_exponent(GS.rescale)
                    .expect("UV momentum power counting failed");
                // The rescaled source already contains the four-dimensional
                // loop measure; retain one energy measure per loop instead.
                let energy_dod = four_d_dod - 3 * loop_count as i32;
                (energy_dod >= 0)
                    .then(|| (cycle.string_label(), loop_count, four_d_dod, energy_dod))
            })
            .collect::<Vec<_>>();
        failures.sort();

        if failures.is_empty() {
            return Ok(());
        }

        let details = failures
            .into_iter()
            .map(|(cycle, loop_count, four_d_dod, energy_dod)| {
                format!("cycle {cycle}: L={loop_count}, DOD_4D={four_d_dod}, DOD_E={energy_dod}")
            })
            .collect::<Vec<_>>()
            .join("; ");
        Err(eyre!(
            "graph `{}` contains a loop-energy nonconvergent integral ({details}); expected DOD_E < 0 for every cycle union",
            self.name,
        ))
    }
}

/// GammaLoop's reference adapter. It keeps the behaviour of the UV code before the
/// move: the requested subgraph is the domain, particles resolve in the caller's
/// model and dummy edges count in particle queries, stored edge and vertex degrees
/// prune, and `compute_dod` classifies with owner-tagged hard momenta.
///
/// The `Result`s are only for the basis query. Particle and power-counting queries
/// never return `Err`; they panic on broken invariants instead, as before the move:
/// particle identifiers must resolve in the caller's model (the graph's canonical
/// model or one sharing its particle table), and the graph numerator must have a
/// valid tensor interface with an integral leading UV power.
impl UvGraph for Graph {
    type Edge = Edge;
    type Vertex = Vertex;
    type HedgeData = HedgeData;
    type Basis = LoopMomentumBasis;
    type Error = LmbError;

    fn hedge_graph(&self) -> &HedgeGraph<Edge, Vertex, HedgeData> {
        &self.underlying
    }

    /// The requested subgraph itself: GammaLoop selects its UV domain (cuts,
    /// initial states and dummy edges) before enumeration.
    fn uv_domain(&self, subgraph: &SuBitGraph) -> SuBitGraph {
        subgraph.clone()
    }

    fn is_dummy_edge(&self, edge: EdgeIndex) -> bool {
        self.underlying[edge].is_dummy
    }

    /// Resolves the particle in the caller's model, dummy edges included.
    ///
    /// # Panics
    ///
    /// Panics if the edge's particle is not in `model`.
    fn edge_pdg_code(&self, edge: EdgeIndex, model: &Model) -> Result<Option<isize>, LmbError> {
        Ok(self.underlying[edge].particle_pdg_code(model))
    }

    /// Whether the edge's particle is massive in the caller's model, dummy edges
    /// included.
    ///
    /// # Panics
    ///
    /// Panics if the edge's particle is not in `model`.
    fn is_massive_edge(&self, edge: EdgeIndex, model: &Model) -> Result<bool, LmbError> {
        Ok(self.underlying[edge].is_massive(model))
    }

    /// Four powers per loop plus the stored edge and vertex degrees.
    fn local_degree_of_divergence(&self, subgraph: &SuBitGraph) -> Result<i32, LmbError> {
        let mut dod: i32 = 4 * self.n_loops(subgraph) as i32;
        for (p, _, e) in self.underlying.iter_edges_of(subgraph) {
            if p.is_paired() {
                dod += e.data.dod.deref();
            }
        }

        for (_, _, n) in self.underlying.iter_nodes_of(subgraph) {
            dod += n.dod.deref();
        }

        Ok(dod)
    }

    /// The numerator-based power counting of [`UltravioletGraph::compute_dod`].
    ///
    /// # Panics
    ///
    /// Panics if the numerator has no valid tensor interface or its leading UV
    /// power is not an integer.
    fn degree_of_divergence(&self, subgraph: &SuBitGraph) -> Result<i32, LmbError> {
        Ok(self.compute_dod(subgraph))
    }

    fn empty_basis(&self) -> LoopMomentumBasis {
        self.empty_lmb()
    }

    fn compatible_sub_basis(
        &self,
        subgraph: &InternalSubGraph,
        parent: &LoopMomentumBasis,
    ) -> Result<LoopMomentumBasis, LmbError> {
        let lmb =
            self.try_compatible_sub_lmb(subgraph, self.dummy_less_full_crown(subgraph), parent)?;
        debug!(
            graph = %self.dot_lmb_of(subgraph, &lmb),
            string_label = %subgraph.string_label(),
            "Routed spinney subgraph"
        );
        Ok(lmb)
    }
}

impl UltravioletGraph for Graph {
    fn denominator<S: SubGraphLike, T: Fn(&Edge) -> isize>(
        &self,
        subgraph: &S,
        model: &Model,
        edge_powers: T,
    ) -> Atom {
        let symbols = feynkit_graph::expressions::PropagatorSymbols {
            momentum: GS.emr_mom,
            denominator: GS.den,
        };
        let result = feynkit_graph::expressions::GraphExpressions::denominator_of(
            &self.underlying,
            subgraph,
            |edge_id, edge| {
                Ok::<_, std::convert::Infallible>(symbols.denominator(
                    edge_id,
                    &edge.mass_atom(model).pow(2),
                    GS.dim,
                ))
            },
            |_, edge| edge_powers(edge),
        );
        match result {
            Ok(denominator) => denominator,
            Err(never) => match never {},
        }
    }

    fn numerator<S: SubGraphLike + SubSetOps>(
        &self,
        subgraph: &S,
        without: &S,
    ) -> Result<Numerator<AppliedFeynmanRule>> {
        let num = Numerator::default();

        num.fill_in_reduced(self, subgraph, without)
    }

    fn compute_dod<S: SubGraphLike<Base = SuBitGraph> + SubSetOps>(&self, subgraph: &S) -> i32 {
        let lmb = self.lmb_of(subgraph);
        let empty = self.underlying.empty_subgraph();
        // Mass terms do not affect the leading UV degree. Building the power-
        // counting denominator with zero masses keeps this structural query
        // independent of model resolution while the physical denominator
        // remains explicitly model-aware.
        let mut denominator = Atom::one();
        for (pair, edge_id, _) in self.underlying.iter_edges_of(subgraph) {
            if pair.is_paired() {
                denominator *= feynkit_graph::expressions::PropagatorSymbols {
                    momentum: GS.emr_mom,
                    denominator: GS.den,
                }
                .denominator(edge_id, &Atom::Zero, GS.dim);
            }
        }
        let integrand = self
            .numerator(subgraph, &empty)
            .expect("graph numerator must have a valid tensor interface")
            .to_d_dim(GS.dim)
            .expect("UV power-counting dimension conversion must preserve tensor interfaces")
            .get_single_atom()
            .unwrap()
            / denominator;
        let nloops: usize = self.n_loops(subgraph);
        self.uv_rescaled(subgraph.included(), nloops, &lmb, &lmb, &integrand)
            .trailing_exponent(GS.rescale)
            .expect("UV momentum power counting failed")
    }
}

pub trait UVE {
    fn mass_atom(&self, model: &crate::model::Model) -> Atom;
    fn particle_pdg_code(&self, model: &crate::model::Model) -> Option<isize>;
    fn is_massive(&self, model: &crate::model::Model) -> bool;
    fn is_fermion(&self, model: &Model) -> bool;
    fn chemical_potential_atom(&self, model: &Model) -> Option<Atom>;
}

#[cfg(test)]
mod shared_expression_tests {
    use super::*;
    use crate::uv::spenso_lor_atom;
    use idenso::tensor::{ContractSettings, SymbolicTensor};
    use linnet::half_edge::involution::EdgeIndex;

    #[test]
    fn shared_propagator_matches_runtime_tensor_contraction() {
        crate::initialisation::test_initialise().unwrap();
        let edge = EdgeIndex(7);
        let mass_squared = symbolica::symbol!("shared_propagator_test::m").pow(2);
        let momentum_square =
            SymbolicTensor::infer_indexed(spenso_lor_atom(edge.0 as i32, edge.0, GS.dim).pow(2))
                .unwrap()
                .contract(ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap()
                .to_dots()
                .unwrap();
        let previous = GS.den(
            edge.0,
            GS.emr_mom.call(edge.0),
            &mass_squared,
            momentum_square.expression() - &mass_squared,
        );
        let shared = feynkit_graph::expressions::PropagatorSymbols {
            momentum: GS.emr_mom,
            denominator: GS.den,
        }
        .denominator(edge, &mass_squared, GS.dim);
        assert_eq!(shared, previous);
    }
}
