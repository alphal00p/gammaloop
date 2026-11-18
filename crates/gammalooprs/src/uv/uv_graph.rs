use std::{cell::RefCell, collections::BTreeSet, ops::Deref};

use ahash::{AHashMap, AHashSet};
use color_eyre::Result;
use eyre::{WrapErr, eyre};
use idenso::color::{ColorSimplifier, ColorSimplifySettings};
use idenso::shorthands::schoonschip::Schoonschip;
use linnet::half_edge::{
    HedgeGraph, PowersetIterator,
    involution::{Flow, Hedge, HedgePair},
    subgraph::{
        Cycle, InternalSubGraph, ModifySubSet, PairwiseSubSetOps, SuBitGraph, SubGraphLike,
        SubGraphOps, SubSetLike, SubSetOps, subset::SubSet,
    },
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    domains::atom::AtomField,
    function,
    poly::series::Series,
};
use thiserror::Error;
use tracing::debug;

use crate::graph::LmbError;
use crate::graph::cuts::CutSet;
use crate::uv::approx::CutStructure;
use crate::uv::thermal_zero_mode::{ZeroModeHostError, assign_zero_mode_hosts};
use crate::{
    graph::{Edge, FeynmanGraph, Graph, HedgeData, LMBext, LoopMomentumBasis, Vertex},
    integrands::process::param_builder::ParamBuilderGraph,
    momentum::sample::LoopIndex,
    numerator::{AppliedFeynmanRule, Numerator, NumeratorStateError},
    settings::global::MediumSettings,
    utils::{GS, W_, symbolica_ext::DOD},
    uv::{ApproximationType, UVgenerationSettings, settings::CTIdentifier},
};

use super::{Spinney, Wood, spenso_lor_atom};

pub trait UltravioletGraph: LMBext + FeynmanGraph + ParamBuilderGraph {
    fn n_loops<S: SubGraphLike, E, V, H>(&self, subgraph: &S) -> usize
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        self.as_ref().cyclotomatic_number(subgraph)
    }

    // fn concretize_spinney<S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> Graph;

    fn dummy_less_full_crown<S: SubGraphLike>(&self, subgraph: &S) -> S::Base
    where
        S::Base: ModifySubSet<Hedge> + SubGraphOps;

    ///Get the numerator of the graph.
    /// If multiply_prefactor is true, the numerator is multiplied by the global  prefactor. (just num not projector)
    fn numerator<S: SubGraphLike + SubSetOps>(
        &self,
        subgraph: &S,
        without: &S,
    ) -> Numerator<AppliedFeynmanRule>;
    fn denominator<S: SubGraphLike, T: Fn(&Edge) -> isize>(
        &self,
        subgraph: &S,
        edge_powers: T,
    ) -> Atom;

    fn boundary_pdg_set<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> BTreeSet<isize>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let graph = self.as_ref();
        graph
            .full_crown(subgraph)
            .included_iter()
            .filter_map(|hedge| {
                let edge_id = graph[&hedge];
                graph[edge_id].particle_pdg_code().map(|pdg| {
                    if graph.flow(hedge) == Flow::Source {
                        -pdg
                    } else {
                        pdg
                    }
                })
            })
            .collect()
    }

    fn internal_pdg_set<E: UVE, V, H, S: SubGraphLike>(&self, subgraph: &S) -> BTreeSet<isize>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        self.as_ref()
            .iter_edges_of(subgraph)
            .filter_map(|(pair, edge_id, _)| {
                pair.is_paired()
                    .then(|| self.as_ref()[edge_id].particle_pdg_code())
                    .flatten()
            })
            .collect()
    }

    fn has_massive_boundary_external<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> bool
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let graph = self.as_ref();
        graph.full_crown(subgraph).included_iter().any(|hedge| {
            let edge_id = graph[&hedge];
            graph[edge_id].is_massive()
        })
    }

    fn ct_identifier<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> CTIdentifier
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        CTIdentifier::new(
            self.boundary_pdg_set(subgraph),
            Some(self.internal_pdg_set(subgraph)),
        )
    }

    fn approximation_scheme<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        settings: &UVgenerationSettings,
        dod: i32,
    ) -> ApproximationType
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        settings.approximation_scheme_for(
            &self.ct_identifier(subgraph),
            dod,
            self.has_massive_boundary_external(subgraph),
        )
    }

    fn classify_spinney<E: UVE, V, H>(
        &self,
        spinney: InternalSubGraph,
        settings: &UVgenerationSettings,
        lmb: &LoopMomentumBasis,
    ) -> Result<Option<Spinney>>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let dod = self.compute_dod(&spinney.filter);
        if dod < 0 {
            return Ok(None);
        }

        let renormalization_scheme = self.approximation_scheme(&spinney.filter, settings, dod);

        if renormalization_scheme != ApproximationType::Unsubtracted {
            Spinney::with_scheme(spinney, self, lmb, renormalization_scheme, dod)
        } else {
            Ok(None)
        }
    }

    fn classified_spinneys<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        settings: &UVgenerationSettings,
        lmb: &LoopMomentumBasis,
    ) -> Result<Vec<Spinney>>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        if !settings.subtract_uv {
            return Ok(vec![Spinney::empty(self)]);
        }

        let graph = self.as_ref();
        let (connected, disconnected): (Vec<_>, Vec<_>) = self
            .spinneys(subgraph)
            .into_iter()
            .partition(|spinney| graph.connected_components(spinney).len() <= 1);

        let mut classified: Vec<_> = connected
            .into_iter()
            .filter_map(|spinney| self.classify_spinney(spinney, settings, lmb).transpose())
            .collect::<Result<_>>()?;
        let retained_factors: AHashSet<_> = classified
            .iter()
            .map(|spinney| spinney.filter().clone())
            .collect();

        for spinney in disconnected {
            let components = graph.connected_components(&spinney);
            if !components.into_iter().all(|component| {
                let factor = InternalSubGraph::cleaned_filter_optimist(component, graph);
                retained_factors.contains(&factor.filter)
            }) {
                continue;
            }

            let dod = self.compute_dod(&spinney.filter);
            if let Some(spinney) =
                Spinney::with_scheme(spinney, self, lmb, ApproximationType::MUV, dod)?
            {
                classified.push(spinney);
            }
        }
        Ok(classified)
    }

    fn all_cycle_unions<E, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> AHashSet<InternalSubGraph>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let ref_graph = self.as_ref();
        let _init_node = ref_graph.iter_nodes_of(subgraph).next().unwrap().0;
        let all_subcycles: Vec<_> =
            Cycle::all_sum_powerset_filter_map(&ref_graph.cycle_basis_of(subgraph).0, &Some)
                .map(|a| a.into_iter().map(|c| c.internal_graph(ref_graph)).collect())
                .unwrap();

        // println!("{}", self.base_dot());
        let spinneys: AHashSet<_> = InternalSubGraph::all_unions_iterative(&all_subcycles);

        spinneys
    }
    fn all_limits<E, V, H, S: SubGraphLike>(
        &self,
        subgraph: &S,
        expr: &Atom,
        expansion: Symbol,
        lmb: &LoopMomentumBasis,
    ) -> Vec<(SubSet<LoopIndex>, Series<AtomField>)>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let mom_reps = self.normal_emr_replacement(subgraph, lmb, &[W_.x___], |_s| true);

        let ose_reps = self.get_ose_replacements();
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

    fn wood<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(&self, subgraph: &S) -> Result<Wood>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        self.wood_with_settings(
            subgraph,
            &UVgenerationSettings::default(),
            &self.as_ref().lmb_of(&self.as_ref().full_filter()),
        )
    }

    fn wood_with_settings<E: UVE, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
        settings: &UVgenerationSettings,
        lmb: &LoopMomentumBasis,
    ) -> Result<Wood>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        Ok(Wood::from_spinneys(
            self.classified_spinneys(subgraph, settings, lmb)?,
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

    fn local_dod<S: SubGraphLike>(&self, subgraph: &S) -> i32;

    /// The static infrared degree `ω(Σ)` of a soft set `Σ`, whose momenta
    /// span `loops` independent loops once the rest of the graph is hard.
    ///
    /// It counts `3` per loop for the static spatial measure, the degree of
    /// every paired edge of `sigma`, and the degree of every vertex all of
    /// whose half-edges are paired and lie in `sigma`, with the edge and
    /// vertex degrees of [`Self::local_dod`]. A vertex touching a hard line
    /// belongs to the hard part. The loop count is an input because a soft
    /// set need not contain its own cycles: the soft lines of a polarization
    /// ring form a path.
    fn static_ir_degree<S: SubGraphLike>(&self, sigma: &S, loops: usize) -> i32;

    fn spinneys<E, V, H, S: SubGraphLike<Base = SuBitGraph>>(
        &self,
        subgraph: &S,
    ) -> AHashSet<InternalSubGraph>
    where
        Self: AsRef<HedgeGraph<E, V, H>>,
    {
        let ref_graph: &HedgeGraph<E, V, H> = self.as_ref();
        debug!(subgraph=%ref_graph.dot(subgraph),"Spinneys of subgraph");
        let _b: SuBitGraph = ref_graph.empty_subgraph();

        if subgraph.is_empty() {
            let mut spinneys = AHashSet::new();
            spinneys.insert(ref_graph.empty_subgraph());
            return spinneys;
        }

        let cycles = ref_graph.cycle_basis_of(subgraph).0;

        let _init_node = ref_graph.iter_nodes_of(subgraph).next().unwrap().0;
        let all_subcycles: Vec<_> = Cycle::all_sum_powerset_filter_map(&cycles, &Some)
            .map(|a| a.into_iter().map(|c| c.internal_graph(ref_graph)).collect())
            .unwrap();

        let dod_cache = RefCell::new(AHashMap::new());

        let mut spinneys: AHashSet<_> = InternalSubGraph::all_ops_iterative_filter_map(
            &all_subcycles,
            &|a, b| a.union(b),
            &|union| {
                // println!("{}", self.as_ref().dot(&union));
                let cached_dod = dod_cache.borrow().get(&union).copied();
                let keep = match cached_dod {
                    Some(keep) => keep,
                    None => {
                        let keep = self.local_dod(&union) >= 0;
                        dod_cache.borrow_mut().insert(union.clone(), keep);
                        keep
                    }
                };
                if keep { Some(union) } else { None }
            },
        );

        spinneys.insert(ref_graph.empty_subgraph());

        spinneys
    }
}

impl AsRef<HedgeGraph<Edge, Vertex, HedgeData>> for Graph {
    fn as_ref(&self) -> &HedgeGraph<Edge, Vertex, HedgeData> {
        &self.underlying
    }
}

impl Graph {
    /// The UV spinneys of every cut of `cuts`, in order: the classified
    /// spinneys of its uncut subgraph, with their thermal zero-mode hosts
    /// assigned and the full observable restricted. Each cut is classified
    /// when the iterator reaches it.
    ///
    /// Zero-mode hosts need a single uncut structure, so restricting each cut
    /// is the restriction of the spinneys of all cuts together.
    ///
    /// # Errors
    ///
    /// Each item fails with the classification, zero-mode host and
    /// full-observable restriction errors of its cut, in the context of the
    /// graph and the cut.
    pub(crate) fn classified_cut_spinneys<'a>(
        &'a self,
        cuts: &'a CutStructure,
        settings: &'a UVgenerationSettings,
        medium: &'a MediumSettings,
    ) -> impl Iterator<Item = Result<Vec<Spinney>>> + 'a {
        cuts.cuts.iter().enumerate().map(move |(cut_index, cut)| {
            self.classified_spinneys_of_cut(cuts, cut_index, cut, settings, medium)
        })
    }

    /// The UV spinneys of `cut`, the cut `cut_index` of `cuts`.
    fn classified_spinneys_of_cut(
        &self,
        cuts: &CutStructure,
        cut_index: usize,
        cut: &CutSet,
        settings: &UVgenerationSettings,
        medium: &MediumSettings,
    ) -> Result<Vec<Spinney>> {
        let cut_subgraph = self
            .full_filter()
            .subtract(&self.initial_state_cut.left)
            .subtract(&cut.union);
        let mut spinneys = self
            .classified_spinneys(&cut_subgraph, settings, &self.loop_momentum_basis)
            .wrap_err_with(|| {
                format!(
                    "graph '{}' cut {cut_index}: failed to classify UV counterterm components",
                    self.name
                )
            })?;
        assign_zero_mode_hosts(self, &cut_subgraph, &mut spinneys, settings, medium, cuts)
            .wrap_err_with(|| {
                format!(
                    "graph '{}' cut {cut_index}: failed to assign thermal zero-mode hosts",
                    self.name
                )
            })?;
        self.restrict_full_observable_spinney(&mut spinneys, medium)
            .wrap_err_with(|| {
                format!(
                    "graph '{}' cut {cut_index}: failed to restrict the full-observable UV counterterm",
                    self.name
                )
            })?;
        Ok(spinneys)
    }

    /// Restrict the operation covering every loop under vacuum subtraction,
    /// and check that every thermal zero-mode host keeps its enclosing
    /// whole-graph U.
    ///
    /// Under vacuum subtraction the full-observable spinney is omitted, unless
    /// some spinney hosts a zero-mode region; it is then kept as a plain U.
    /// Without vacuum subtraction the classification is unchanged.
    ///
    /// # Errors
    ///
    /// In either vacuum-subtraction mode, a zero-mode host whose full
    /// observable is not retained: the host's static cograph is
    /// log-divergent at its top Taylor order, and only the whole-graph U
    /// subtracts it.
    pub(crate) fn restrict_full_observable_spinney(
        &self,
        spinneys: &mut Vec<Spinney>,
        medium: &MediumSettings,
    ) -> Result<(), ZeroModeHostError> {
        // Spinneys contain no bridges.
        let full_observable = InternalSubGraph::cleaned_filter_pessimist(
            self.underlying.non_bridges(),
            self.as_ref(),
        );
        let has_zero_mode_host = spinneys
            .iter()
            .any(|spinney| spinney.thermal_zero_mode().is_some());
        if medium.vacuum_subtraction {
            // Completed thermal weights already carry (1 - V). An operation
            // covering all loops leaves no medium dependence, even when tree
            // edges survive, so its U vanishes on every branch that projects
            // all occupations. A frozen zero-mode branch is a vertex for every
            // enclosing operation, so the whole-graph U is the UV subtraction
            // of its static cograph. Keep that U alone: without external
            // momenta, the whole graph's soft expansion would reproduce the
            // integrand.
            if has_zero_mode_host {
                // The full observable leaves no soft set, so it hosts nothing.
                for spinney in spinneys
                    .iter_mut()
                    .filter(|spinney| spinney.subgraph == full_observable)
                {
                    spinney.demote_to_muv();
                }
            } else {
                // Omit it before either 4D or 3D computation.
                spinneys.retain(|spinney| {
                    spinney.subgraph.is_empty() || spinney.subgraph != full_observable
                });
            }
        }

        if let Some(host) = spinneys
            .iter()
            .find(|spinney| spinney.thermal_zero_mode().is_some())
            && !spinneys
                .iter()
                .any(|spinney| spinney.subgraph == full_observable)
        {
            return Err(ZeroModeHostError::MissingEnclosingSpinney {
                host: self
                    .underlying
                    .paired_edges(host.filter())
                    .into_iter()
                    .collect(),
                dod: self.compute_dod(&full_observable.filter),
            });
        }
        Ok(())
    }

    /// The graph-canonical loop-momentum basis of `subgraph`: the route the
    /// graph's own basis induces on it, with the crown of the subgraph less
    /// dummy edges as its external flows. The whole internal graph keeps the
    /// graph's basis. Every Taylor chart of a direct local-3D operation, and
    /// the coordinates of a thermal zero-mode host, start from this basis.
    ///
    /// # Errors
    ///
    /// The graph's basis induces no basis on `subgraph`.
    pub(crate) fn canonical_component_lmb(
        &self,
        subgraph: &SuBitGraph,
    ) -> Result<LoopMomentumBasis, LmbError> {
        let external: SuBitGraph = self.external_filter();
        let full_internal = self.full_filter().subtract(&external);
        if *subgraph == full_internal
            && self.n_loops(&full_internal) == self.loop_momentum_basis.loop_edges.len()
        {
            return Ok(self.loop_momentum_basis.clone());
        }
        self.try_compatible_sub_lmb(
            subgraph,
            self.dummy_less_full_crown(subgraph),
            &self.loop_momentum_basis,
        )
    }

    /// The Feynman-rule numerator that the edges of `subgraph` add to
    /// `given`, with its colour algebra closed before a residue map repeats
    /// it across rows. Its Lorentz and spinor structure stays as it is.
    ///
    /// # Errors
    ///
    /// The applied numerator is not a single expression.
    pub(crate) fn color_reduced_numerator(
        &self,
        subgraph: &SuBitGraph,
        given: &SuBitGraph,
    ) -> Result<Atom, NumeratorStateError> {
        Ok(self
            .numerator(subgraph, given)
            .get_single_atom()?
            .simplify_color_with(ColorSimplifySettings {
                simplify_non_color: false,
                ..Default::default()
            }))
    }

    /// The degree of the measures of `loops` loops, `per_loop` each.
    fn loop_measure_degree(per_loop: i32, loops: usize) -> i32 {
        per_loop * i32::try_from(loops).expect("a graph has fewer than 2^31 loops")
    }

    /// The degree of the paired edges of `subgraph`: its propagators and
    /// their numerators.
    fn paired_edge_degree<S: SubGraphLike>(&self, subgraph: &S) -> i32 {
        self.underlying
            .iter_edges_of(subgraph)
            .filter(|(pair, _, _)| pair.is_paired())
            .map(|(_, _, edge)| *edge.data.dod.deref())
            .sum()
    }

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
                let integrand = &numerator / self.denominator(&cycle.filter, |_| 1);
                let four_d_dod = self
                    .uv_rescaled(cycle.filter.included(), loop_count, &lmb, &lmb, &integrand)
                    .trailing_exponent();
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

impl UltravioletGraph for Graph {
    fn dummy_less_full_crown<S: SubGraphLike>(&self, subgraph: &S) -> S::Base
    where
        S::Base: ModifySubSet<Hedge>,
    {
        let a = self.full_crown(subgraph);
        let mut ac = a.clone();

        a.included_iter().for_each(|a| {
            if self[self[&a]].is_dummy {
                ac.sub(a);
            }
        });

        ac
    }

    // fn concretize_spinney<S: SubGraphLike<Base = SuBitGraph>>(&self, s: &S) -> Graph {
    //     let mut underlying = self.underlying.clone();
    //     let dod = self.dod(s);
    //     let cc = underlying.connected_components(s);

    //     for c in cc{

    //     let v = Vertex {

    //     }
    //     underlying.identify_nodes_without_self_edges(s, node_data_merge);
    //     };
    //     }

    //     self.clone()
    // }

    fn denominator<S: SubGraphLike, T: Fn(&Edge) -> isize>(
        &self,
        subgraph: &S,
        edge_powers: T,
    ) -> Atom {
        let mut den = Atom::num(1);

        for (pair, eid, d) in self.underlying.iter_edges_of(subgraph) {
            if matches!(pair, HedgePair::Paired { .. }) {
                let m2 = d.data.mass_atom().pow(2);
                let edge_power = edge_powers(d.data);
                let is_power_negative = edge_power < 0;
                let prop_den = GS.den(
                    usize::from(eid),
                    function!(GS.emr_mom, usize::from(eid)),
                    &m2,
                    spenso_lor_atom(usize::from(eid) as i32, usize::from(eid), GS.dim)
                        .pow(2)
                        .to_dots()
                        - &m2,
                );
                for _i in 0..edge_power.abs() {
                    if is_power_negative {
                        den /= prop_den.clone();
                    } else {
                        den *= prop_den.clone();
                    }
                }
            }
        }

        den
    }
    fn numerator<S: SubGraphLike + SubSetOps>(
        &self,
        subgraph: &S,
        without: &S,
    ) -> Numerator<AppliedFeynmanRule> {
        let num = Numerator::default();

        num.fill_in_reduced(self, subgraph, without)
    }

    fn compute_dod<S: SubGraphLike<Base = SuBitGraph> + SubSetOps>(&self, subgraph: &S) -> i32 {
        let lmb = self.lmb_of(subgraph);
        let empty = self.underlying.empty_subgraph();
        let integrand = self
            .numerator(subgraph, &empty)
            .to_d_dim(GS.dim)
            .get_single_atom()
            .unwrap()
            / self.denominator(subgraph, |_| 1);
        let nloops: usize = self.n_loops(subgraph);
        self.uv_rescaled(subgraph.included(), nloops, &lmb, &lmb, &integrand)
            .trailing_exponent()
    }

    fn local_dod<S: SubGraphLike>(&self, subgraph: &S) -> i32 {
        let mut dod = Self::loop_measure_degree(4, self.n_loops(subgraph))
            + self.paired_edge_degree(subgraph);
        for (_, _, n) in self.underlying.iter_nodes_of(subgraph) {
            dod += n.dod.deref();
        }

        dod
    }

    fn static_ir_degree<S: SubGraphLike>(&self, sigma: &S, loops: usize) -> i32 {
        let mut degree = Self::loop_measure_degree(3, loops) + self.paired_edge_degree(sigma);
        for (_, mut crown, vertex) in self.underlying.iter_nodes_of(sigma) {
            if crown.all(|hedge| {
                let partner = self.underlying.inv(hedge);
                partner != hedge && sigma.includes(&hedge) && sigma.includes(&partner)
            }) {
                degree += vertex.dod.deref();
            }
        }

        degree
    }
}

/// Why a line has no Matsubara zero mode.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Error)]
pub enum NoZeroMode {
    #[error("is a fermion, whose antiperiodic frequencies have no zero mode")]
    Fermion,
    #[error("has a chemical potential, which shifts its frequencies off zero")]
    ChemicalPotential,
}

pub trait UVE {
    fn mass_atom(&self) -> Atom;
    fn particle_pdg_code(&self) -> Option<isize>;
    fn is_massive(&self) -> bool;
    fn is_fermion(&self) -> bool;
    fn chemical_potential_atom(&self) -> Option<Atom>;

    /// Whether the line has a Matsubara zero mode: only a boson without a
    /// chemical potential has periodic, unshifted frequencies, one of which
    /// vanishes.
    ///
    /// # Errors
    ///
    /// Why the line has none.
    fn zero_mode(&self) -> Result<(), NoZeroMode> {
        if self.is_fermion() {
            Err(NoZeroMode::Fermion)
        } else if self.chemical_potential_atom().is_some() {
            Err(NoZeroMode::ChemicalPotential)
        } else {
            Ok(())
        }
    }
}
