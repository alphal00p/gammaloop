use std::{
    collections::HashSet,
    fmt,
    sync::{Arc, OnceLock},
};

use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use colored::Colorize;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionArgument, FunctionBuilder, SliceType, Symbol},
    function,
    id::{Condition, Match},
    solve::SolveCoverage,
};

use crate::{
    EvaluationOrder, Integral, ReplacementRules, VakintError, function_condition,
    get_individual_momenta, get_individual_momenta_from_atom, get_prop_with_id, graph::Graph,
    symbols::S,
};
mod catalogue;
pub(crate) mod five_loop;
pub(crate) use catalogue::FiveLoopCatalogue;
pub use catalogue::FiveLoopCatalogueStatistics;
mod momentum_mapping;
pub use momentum_mapping::MomentumMatchLimits;

pub enum TopologyContractions {
    Custom(Vec<Vec<usize>>),
    Automatic,
}

#[derive(Debug, Clone)]
pub struct Topologies(Vec<Topology>, OnceLock<Arc<FiveLoopCatalogue>>, Arc<()>);

/// A retained incarnation, not a graph hash or an allocation-address integer.
/// Clones share the immutable base registry; replacing the public registry
/// gives a different identity. The optional catalogue is installed separately.
#[derive(Debug, Clone)]
pub(crate) struct RegistryIdentity {
    base: Arc<()>,
    catalogue: Option<Arc<FiveLoopCatalogue>>,
}

impl PartialEq for RegistryIdentity {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.base, &other.base)
            && match (&self.catalogue, &other.catalogue) {
                (Some(left), Some(right)) => Arc::ptr_eq(left, right),
                (None, None) => true,
                _ => false,
            }
    }
}

impl Eq for RegistryIdentity {}

impl fmt::Display for Topologies {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            self.iter()
                .enumerate()
                .map(|(i, t)| format!("{}) {}", i + 1, t))
                .collect::<Vec<_>>()
                .join("\n")
        )
    }
}

impl Topologies {
    pub(crate) fn registry_identity(&self) -> RegistryIdentity {
        RegistryIdentity {
            base: Arc::clone(&self.2),
            catalogue: self.1.get().cloned(),
        }
    }

    /// Construct the ordinary one-through-four-loop registry. The saved
    /// five-loop catalogue is installed only by successful artifact loading.
    pub fn generate_topologies() -> Result<Self, VakintError> {
        // ==
        // One-loop topology
        // ==
        let mut topologies = Topologies(
            vec![
                Integral::new(
                    1,
                    Some(vk_parse!("topo(prop(1,edge(1,1),k(1),msq(1),pow(1)))").unwrap()),
                    Some(vk_parse!("I1L(msq(1),pow(1))").unwrap()),
                    EvaluationOrder::all_but_fmft(),
                )?
                .into(),
            ],
            OnceLock::new(),
            Arc::new(()),
        );
        // ==
        // Two-loop topologies
        // ==
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                vk_parse!(
                    "topo(\
                         prop(1,edge(1,2),k(1),msq(1),pow(1))\
                        *prop(2,edge(1,2),k(2),msq(1),pow(2))\
                        *prop(3,edge(2,1),k(1)+k(2),msq(1),pow(3))\
                    )"
                )
                .unwrap()
                .as_view(),
                vk_parse!("I2L(msq(1),pow(1),pow(2),pow(3))")
                    .unwrap()
                    .as_view(),
                TopologyContractions::Custom(vec![vec![], vec![3]]),
                EvaluationOrder::all_but_fmft(),
            )?
            .0,
        );
        // ==
        // Three-loop topologies
        // ==
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                vk_parse!(
                    "topo(\
                         prop(1,edge(1,2),k(1),msq(1),pow(1))\
                        *prop(2,edge(2,3),k(2),msq(1),pow(2))\
                        *prop(3,edge(3,1),k(3),msq(1),pow(3))\
                        *prop(4,edge(1,4),k(3)-k(1),msq(1),pow(4))\
                        *prop(5,edge(2,4),k(1)-k(2),msq(1),pow(5))\
                        *prop(6,edge(3,4),k(2)-k(3),msq(1),pow(6))\
                    )"
                )
                .unwrap()
                .as_view(),
                vk_parse!("I3L(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6))")
                    .unwrap()
                    .as_view(),
                // These contraction IDs are established public short forms.
                // Keep their representatives stable when graph quotient or
                // canonical-labeling implementations choose another isomorph.
                TopologyContractions::Custom(vec![
                    vec![],
                    vec![6],
                    vec![1, 6],
                    vec![3, 6],
                    vec![1, 3, 6],
                ]),
                EvaluationOrder::all_but_fmft(),
            )?
            .0,
        );
        // ==
        // Four-loop topologies
        // ==
        // H topology in FMFT
        // in FMFT paper, FMFT edge ids are these edge ids + 1, except for p_1 which is edge #9 here.
        // (FMFT -> Vakint) : (1->9, 2->1, 3->2, 4->3, 5->4, 6->5, 7->6, 8->7, 9->8)
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                vk_parse!(
                    "topo(\
                         prop(1,edge(5,1),k(1),msq(1),pow(1))\
                        *prop(2,edge(2,6),k(2),msq(1),pow(2))\
                        *prop(3,edge(6,5),k(3),msq(1),pow(3))\
                        *prop(4,edge(3,4),k(4),msq(1),pow(4))\
                        *prop(5,edge(4,5),k(1)-k(3),msq(1),pow(5))\
                        *prop(6,edge(6,3),k(2)-k(3),msq(1),pow(6))\
                        *prop(7,edge(4,1),k(3)-k(1)+k(4),msq(1),pow(7))\
                        *prop(8,edge(2,3),k(3)-k(2)+k(4),msq(1),pow(8))\
                        *prop(9,edge(1,2),k(3)+k(4),msq(1),pow(9))\
                    )"
                )
                .unwrap()
                .as_view(),
                vk_parse!(
                    "I4L_H(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6),pow(7),pow(8),pow(9))"
                )
                .unwrap()
                .as_view(),
                TopologyContractions::Custom(vec![vec![]]),
                EvaluationOrder::fmft_and_numerical(None, None),
            )?
            .0,
        );
        // X topology in FMFT
        // in FMFT paper, FMFT edge ids are these edge ids + 1.
        // (FMFT -> Vakint) : (2->1, 3->2, 4->3, 5->4, 6->5, 7->6, 8->7, 9->8, 10->9)
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                vk_parse!(
                    "topo(\
                         prop(1,edge(5,1),k(1),msq(1),pow(1))\
                        *prop(2,edge(2,6),k(2),msq(1),pow(2))\
                        *prop(3,edge(6,5),k(3),msq(1),pow(3))\
                        *prop(4,edge(4,3),k(4),msq(1),pow(4))\
                        *prop(5,edge(3,5),k(1)-k(3),msq(1),pow(5))\
                        *prop(6,edge(6,4),k(2)-k(3),msq(1),pow(6))\
                        *prop(7,edge(3,2),k(3)-k(1)+k(4),msq(1),pow(7))\
                        *prop(8,edge(1,4),k(3)-k(2)+k(4),msq(1),pow(8))\
                        *prop(9,edge(2,1),k(3)-k(1)-k(2)+k(4),msq(1),pow(9))\
                    )"
                )
                .unwrap()
                .as_view(),
                vk_parse!(
                    "I4L_X(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6),pow(7),pow(8),pow(9))"
                )
                .unwrap()
                .as_view(),
                TopologyContractions::Custom(vec![vec![]]),
                EvaluationOrder::fmft_and_numerical(None, None),
            )?
            .0,
        );
        // BMW topology in FMFT
        // When sorting both these edges and the FMT ones according to their ID, both lists match.
        // (FMFT -> Vakint) : (3->1, 4->2, 5->3, 6->4, 7->5, 8->6, 9->7, 10->8)
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                vk_parse!(
                    "topo(\
                         prop(1,edge(1,2),k(1),msq(1),pow(1))\
                        *prop(2,edge(2,5),k(2),msq(1),pow(2))\
                        *prop(3,edge(3,4),k(3),msq(1),pow(3))\
                        *prop(4,edge(4,5),k(4),msq(1),pow(4))\
                        *prop(5,edge(2,3),k(1)-k(2),msq(1),pow(5))\
                        *prop(6,edge(4,1),k(3)-k(4),msq(1),pow(6))\
                        *prop(7,edge(5,3),k(2)+k(3)-k(1),msq(1),pow(7))\
                        *prop(8,edge(1,5),k(3)-k(4)-k(1),msq(1),pow(8))\
                    )"
                )
                .unwrap()
                .as_view(),
                vk_parse!(
                    "I4L_BMW(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6),pow(7),pow(8))"
                )
                .unwrap()
                .as_view(),
                TopologyContractions::Custom(vec![vec![]]),
                EvaluationOrder::fmft_and_numerical(None, None),
            )?
            .0,
        );
        // FG topology in FMFT
        // When sorting both these edges and the FMT ones according to their ID, both lists match.
        // (FMFT -> Vakint) : (1->1, 3->2, 4->3, 5->4, 6->5, 7->6, 8->7, 9->8)
        topologies.0.extend(
            Topologies::generate_topologies_with_contractions(
                vk_parse!(
                    "topo(\
                         prop(1,edge(5,3),k(1),msq(1),pow(1))\
                        *prop(2,edge(3,4),k(2),msq(1),pow(2))\
                        *prop(3,edge(4,5),k(3),msq(1),pow(3))\
                        *prop(4,edge(2,1),k(1)-k(3),msq(1),pow(4))\
                        *prop(5,edge(5,1),k(4),msq(1),pow(5))\
                        *prop(6,edge(4,2),k(2)-k(3),msq(1),pow(6))\
                        *prop(7,edge(1,5),k(1)-k(3)+k(4),msq(1),pow(7))\
                        *prop(8,edge(3,2),k(1)-k(2),msq(1),pow(8))\
                    )"
                )
                .unwrap()
                .as_view(),
                vk_parse!("I4L_FG(msq(1),pow(1),pow(2),pow(3),pow(4),pow(5),pow(6),pow(7),pow(8))")
                    .unwrap()
                    .as_view(),
                // Place all further pinched topologies on that FX one
                TopologyContractions::Automatic, //TopologyContractions::Custom(vec![vec![], vec![3]]),
                EvaluationOrder::fmft_and_numerical(None, None),
            )?
            .0,
        );

        topologies.0.push(Topology::Unknown(Integral::new(
            0,
            None,
            None,
            EvaluationOrder::pysecdec_only(None),
        )?));

        Ok(topologies)
    }

    pub fn count_propagators_in_integral(integral: AtomView) -> usize {
        let mut n_props: usize = 0;
        let p = vk_parse!("prop(args__)").unwrap().to_pattern();
        let rhs_p = S.one.to_pattern();
        let matcher = integral.replace(&p);
        let mut match_iter = matcher.iter(&rhs_p);
        while match_iter.next().is_some() {
            n_props += 1;
        }
        n_props
    }

    pub fn generate_topologies_with_contractions(
        canonical_expression: AtomView,
        short_expression: AtomView,
        contractions: TopologyContractions,
        applicable_evaluation_methods: EvaluationOrder,
    ) -> Result<Self, VakintError> {
        let n_tot_props = Topologies::count_propagators_in_integral(canonical_expression);
        if n_tot_props == 0 {
            return Err(VakintError::InvalidIntegralFormat(format!(
                "No propagators found in integral '{}'",
                canonical_expression
            )));
        }
        let mut topologies = vec![];

        let unique_contractions = match contractions {
            TopologyContractions::Custom(c) => c,
            TopologyContractions::Automatic => {
                let master_topology = Topology::generate_topology_with_contraction(
                    n_tot_props,
                    canonical_expression,
                    short_expression,
                    vec![],
                    applicable_evaluation_methods.clone(),
                )?;
                let cs = master_topology
                    .get_integral()
                    .graph
                    .find_unique_contractions();
                topologies.push(master_topology);
                cs
            }
        };

        for contracted_prop_indices in unique_contractions {
            topologies.push(Topology::generate_topology_with_contraction(
                n_tot_props,
                canonical_expression,
                short_expression,
                contracted_prop_indices,
                applicable_evaluation_methods.clone(),
            )?);
        }
        Ok(Topologies(topologies, OnceLock::new(), Arc::new(())))
    }

    /// Match the currently available base/custom and explicitly loaded data.
    /// This registry-only API never prepares a catalogue or opens an artifact.
    pub fn match_topologies_to_user_input(
        &self,
        input: AtomView,
        allow_unknown_integrals: bool,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        let match_started = std::time::Instant::now();
        // The optional five-loop parser must not narrow the legacy Unknown
        // contract before any five-loop matching data has been installed.
        if self
            .iter()
            .any(|topology| matches!(topology, Topology::FiveLoop(_)))
            && let Some(source) = momentum_mapping::MomentumInput::from_expression(input, 5)?
        {
            log::debug!(target: "vakint::five_loop_match", "phase=native_input_ready elapsed_seconds={:.6}", match_started.elapsed().as_secs_f64());
            if let Some(replacement) = source.match_topologies(
                self.iter()
                    .filter(|topology| matches!(topology, Topology::FiveLoop(_))),
            )? {
                return Ok(Some(replacement));
            }
        } else {
            // A validated five-loop long form cannot be a lower-loop topology.
            // In particular factorized five-loop incidence can otherwise trigger
            // expensive unsuccessful lower-loop pattern permutations.
            log::debug!(target: "vakint::five_loop_match", "phase=legacy_sweep_begin");
            for topology in self.iter() {
                if matches!(topology, Topology::Unknown(_)) {
                    continue;
                }
                if let Topology::FiveLoop(integral) = topology {
                    if let Some(mut replacement) =
                        integral.match_integral_to_short_user_input(input)?
                    {
                        replacement.canonical_topology = topology.clone();
                        return Ok(Some(replacement));
                    }
                    // Five-loop integral classes can have different graph incidence.
                    // Their long form has one native verified momentum-mapping path.
                    continue;
                }
                if let Some(replacement_rules) = topology.match_topology_to_user_input(input)? {
                    return Ok(Some(replacement_rules));
                }
            }
            log::debug!(target: "vakint::five_loop_match", "phase=legacy_sweep_end elapsed_seconds={:.6}", match_started.elapsed().as_secs_f64());
        }
        if allow_unknown_integrals {
            for topology in self
                .0
                .iter()
                .filter(|topology| matches!(topology, Topology::Unknown(_)))
            {
                if let Some(replacement) = topology.match_topology_to_user_input(input)? {
                    return Ok(Some(replacement));
                }
            }
        }
        Ok(None)
    }

    /// Scope-aware exact witness selection for long-form five-loop inputs.
    /// Ordinary and lower-loop matching retain their existing fast behavior.
    /// `None` from the criterion means only an ordinary scope rejection; typed
    /// errors propagate. If no witness is admitted, retain an exact witness so
    /// the caller can attempt a partial reduction outside recorded coverage.
    pub fn match_topologies_to_user_input_with_acceptance<K: Ord>(
        &self,
        input: AtomView,
        allow_unknown_integrals: bool,
        limits: MomentumMatchLimits,
        mut criterion: impl FnMut(&ReplacementRules) -> Result<Option<K>, VakintError>,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        if self
            .iter()
            .any(|topology| matches!(topology, Topology::FiveLoop(_)))
            && let Some(source) = momentum_mapping::MomentumInput::from_expression(input, 5)?
            && let Some(candidate) = source.match_topologies_with_acceptance(
                self.iter()
                    .filter(|topology| matches!(topology, Topology::FiveLoop(_))),
                limits,
                &mut criterion,
            )?
        {
            return Ok(Some(candidate));
        }
        let matched = self.match_topologies_to_user_input(input, allow_unknown_integrals)?;
        if let Some(candidate) = &matched
            && matches!(candidate.canonical_topology, Topology::FiveLoop(_))
        {
            // The unchanged parametric/short-form fallback may supply one
            // extra exact witness, but never suppress a criterion error.
            let _ = criterion(candidate)?;
        }
        Ok(matched)
    }
}

#[allow(unused)]
#[derive(Debug, Clone)]
pub enum Topology {
    OneLoop(Integral),
    TwoLoop(Integral),
    ThreeLoop(Integral),
    FourLoop(Integral),
    FiveLoop(Integral),
    Unknown(Integral),
}

impl fmt::Display for Topology {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Topology::OneLoop(i) => write!(f, "{}{}{}", "OneLoop(".magenta(), i, ")".magenta()),
            Topology::TwoLoop(i) => write!(f, "{}{}{}", "TwoLoop(".magenta(), i, ")".magenta()),
            Topology::ThreeLoop(i) => write!(f, "{}{}{}", "ThreeLoop(".magenta(), i, ")".magenta()),
            Topology::FourLoop(i) => write!(f, "{}{}{}", "FourLoop(".magenta(), i, ")".magenta()),
            Topology::FiveLoop(i) => write!(f, "{}{}{}", "FiveLoop(".magenta(), i, ")".magenta()),
            Topology::Unknown(i) => write!(f, "{}{}{}", "Unknown(".red(), i, ")".red()),
        }
    }
}

impl From<Integral> for Topology {
    fn from(i: Integral) -> Self {
        match i.n_loops {
            1 => Topology::OneLoop(i),
            2 => Topology::TwoLoop(i),
            3 => Topology::ThreeLoop(i),
            4 => Topology::FourLoop(i),
            5 => Topology::FiveLoop(i),
            _ => Topology::Unknown(i),
        }
    }
}

#[cfg(test)]
#[path = "topologies/routing_tests.rs"]
mod routing_tests;

impl Topology {
    /// The same exact head is used by native contraction generation and the
    /// cheap optional-catalogue gate; input names are never prefix-matched.
    pub(super) fn contracted_short_symbol(head: Symbol, contraction: &[usize]) -> Symbol {
        if contraction.is_empty() {
            head
        } else {
            vk_symbol!(format!(
                "{}_pinch_{}",
                head,
                contraction
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join("_")
            ))
        }
    }

    pub fn generate_topology_with_contraction(
        n_tot_props: usize,
        canonical_expression: AtomView,
        short_expression: AtomView,
        contraction: Vec<usize>,
        applicable_evaluation_methods: EvaluationOrder,
    ) -> Result<Topology, VakintError> {
        let mut contracted_canonical_expression = canonical_expression.to_owned();
        let mut contracted_short_expression = short_expression.to_owned();

        for prop_id in contraction.iter() {
            if get_prop_with_id(contracted_canonical_expression.as_view(), *prop_id).is_none() {
                return Err(VakintError::InvalidIntegralFormat(format!(
                    "Cannot contract propagatore id {} as it is not found in integral: {}",
                    prop_id, contracted_canonical_expression
                )));
            }
            contracted_canonical_expression = contracted_canonical_expression
                .replace(
                    vk_parse!(format!("prop({},args__)", prop_id).as_str())
                        .unwrap()
                        .to_pattern(),
                )
                .with(S.one.to_pattern());

            for a in [&mut contracted_short_expression] {
                if !a.is_zero() {
                    *a = a
                        .replace(
                            vk_parse!(format!("pow({})", prop_id).as_str())
                                .unwrap()
                                .to_pattern(),
                        )
                        .with(S.zero.to_pattern());
                }
            }
        }

        let mut contracted_graph = Graph::new_from_atom(canonical_expression, n_tot_props)?;
        contracted_graph.contract_edges(&contraction.iter().copied().collect());
        for edge in contracted_graph.edges.values() {
            contracted_canonical_expression = contracted_canonical_expression
                .replace(
                    vk_parse!(format!("prop({},edge(nl_,nr_),args__)", edge.id))
                        .unwrap()
                        .to_pattern(),
                )
                .with(
                    vk_parse!(format!(
                        "prop({},edge({},{}),args__)",
                        edge.id, edge.left_node_id, edge.right_node_id
                    ))
                    .unwrap()
                    .to_pattern(),
                );
        }

        // let tt: Topology = Integral::new(
        //     n_tot_props,
        //     Some(contracted_canonical_expression.clone()),
        //     Some(contracted_short_expression.clone()),
        //     applicable_evaluation_methods.clone(),
        // )?
        // .into();
        // println!("tt={}", tt);
        let (rotated_expression, parent_coordinates) =
            Self::canonicalize_lmb(contracted_canonical_expression.as_view(), n_tot_props)?;
        contracted_canonical_expression = rotated_expression;
        // if contracted_short_expression.to_canonical_string()
        //     == String::from("I3L(msq(1),0,pow(2),pow(3),pow(4),pow(5),0)")
        // {
        //     println!("after={}", tt);
        //     panic!("TTTT");
        // }
        if !contraction.is_empty() {
            contracted_short_expression = contracted_short_expression
                .replace(vk_parse!("TopoName_(args___)").unwrap().to_pattern())
                .when(Condition::from((
                    vk_symbol!("TopoName_"),
                    function_condition(),
                )))
                .with_map(move |match_in| {
                    let topo_name = if let Match::FunctionName(f) =
                        match_in.get(vk_symbol!("TopoName_")).unwrap()
                    {
                        *f
                    } else {
                        panic!("Topology name must correspond to a function name.");
                    };
                    let args = if let Match::Multiple(SliceType::Arg, aviews) =
                        match_in.get(vk_symbol!("args___")).unwrap()
                    {
                        aviews
                    } else {
                        panic!("Topology arguments must be wrapped in Symbolica Arg construction.");
                    };
                    let mut res = FunctionBuilder::new(Self::contracted_short_symbol(
                        topo_name,
                        &contraction,
                    ));
                    for a in args.iter() {
                        res = a.to_owned().add_arg_to_function_builder(res);
                    }
                    res.finish()
                });
        }
        let mut integral = Integral::new(
            n_tot_props,
            Some(contracted_canonical_expression),
            Some(contracted_short_expression),
            applicable_evaluation_methods,
        )?;
        integral.parent_routing = Some((canonical_expression.to_owned(), parent_coordinates));
        Ok(integral.into())
    }

    fn solve_loop_coordinates(
        system: &[Atom],
        variables: &[Atom],
    ) -> Result<Vec<Atom>, VakintError> {
        // The general solve API also accepts nonlinear equations. Retain the
        // native linearity gate before asking it for a complete routing branch.
        Atom::system_to_matrix::<u8, _, _>(system, variables).map_err(|error| {
            VakintError::InvalidIntegralFormat(format!(
                "Could not form the linear system to force the loop momentum basis: {error:?}"
            ))
        })?;
        let solutions = Atom::solve(system)
            .wrt_with_exponent::<u8, _>(variables)
            .map_err(|error| {
                VakintError::InvalidIntegralFormat(format!(
                    "Could not solve the linear system to force the loop momentum basis: {error:?}"
                ))
            })?;
        // A routing witness must be unique and unconditional. Free coordinates
        // or exceptional parameter conditions cannot be dropped when changing
        // the numerator and every propagator simultaneously.
        if solutions.coverage() != SolveCoverage::Complete || !solutions.coverage_guard().is_empty()
        {
            return Err(VakintError::InvalidIntegralFormat(format!(
                "Loop momentum basis requires complete guard-free coverage: {solutions:?}"
            )));
        }
        if solutions.len() != 1 {
            return Err(VakintError::InvalidIntegralFormat(format!(
                "Expected one loop momentum basis, found {} solution branches",
                solutions.len()
            )));
        }
        let solution = &solutions[0];
        if !solution.is_point() || solution.coordinates().len() != variables.len() {
            return Err(VakintError::InvalidIntegralFormat(format!(
                "Loop momentum basis is not unique and unconditional: {solution:?}"
            )));
        }
        Ok(solution
            .coordinates()
            .iter()
            .map(|(_, value)| value.clone())
            .collect())
    }

    pub fn force_an_lmb(
        input_canonical_expression: AtomView,
        n_props: usize,
    ) -> Result<Atom, VakintError> {
        Self::canonicalize_lmb(input_canonical_expression, n_props)
            .map(|(expression, _)| expression)
    }

    fn canonicalize_lmb(
        input_canonical_expression: AtomView,
        n_props: usize,
    ) -> Result<(Atom, Arc<[Atom]>), VakintError> {
        let mut loop_mom_ids: HashSet<i64> = HashSet::new();
        let mut momenta = vec![];
        let mut prop_ids = vec![];
        for prop_id in 1..=n_props {
            if let Some(m) = get_prop_with_id(input_canonical_expression, prop_id) {
                prop_ids.push(prop_id);
                let q = m.get(&vk_symbol!("q_")).unwrap().to_owned();
                for (_s, (_a_id, mom_id)) in get_individual_momenta_from_atom(q.as_view())?.iter() {
                    loop_mom_ids.insert(*mom_id);
                }
                momenta.push((prop_id, q));
            }
        }
        if loop_mom_ids != HashSet::from_iter(1..=loop_mom_ids.len() as i64) {
            return Err(VakintError::InvalidIntegralFormat(format!(
                "Loop momenta must have IDs ranging from 1 to the number of loop momenta: {:?}",
                loop_mom_ids
            )));
        }
        let n_loops = loop_mom_ids.len();

        let mut need_to_force_lmb = false;
        for i_loop in 1..=n_loops {
            if !momenta.iter().any(|(_i_prop, k)| {
                k.as_view() == function!(S.k, Atom::num(i_loop as i64)).as_view()
            }) {
                need_to_force_lmb = true;
                break;
            }
        }
        // if !need_to_force_lmb && False {
        if !need_to_force_lmb {
            let coordinates = (1..=n_loops)
                .map(|index| function!(S.k, Atom::num(index)))
                .collect::<Vec<_>>();
            return Ok((input_canonical_expression.to_owned(), coordinates.into()));
        }

        let g = Graph::new_from_atom(input_canonical_expression, n_props)?;

        let lmb = g.get_one_lmb()?;
        // println!("input_canonical_expression={}", input_canonical_expression);
        // println!("lmb: {:?}", lmb);
        if lmb.len() != n_loops {
            return Err(VakintError::InvalidIntegralFormat(format!(
                "Loop momentum basis identified does not have the same loop count ({}) as detected for this topology ({}): {}.",
                lmb.len(),
                n_loops,
                g
            )));
        }

        let mom_vecs_to_symbols = Arc::new(
            (1..=n_loops)
                .map(|i| {
                    (
                        function!(S.k, Atom::num(i as i64)).to_pattern(),
                        (
                            Atom::var(vk_symbol!(format!("k{}", i))).to_pattern(),
                            (
                                vk_symbol!(format!("k{}", i)),
                                vk_symbol!(format!("krotated{}", i)),
                            ),
                        ),
                    )
                })
                .collect::<Vec<_>>(),
        );

        let mut system = vec![];
        let mut parent_coordinates = Vec::with_capacity(n_loops);
        for (i_lmb, lmb_edge_id) in lmb.iter().enumerate() {
            if let Some(m) = get_prop_with_id(input_canonical_expression, *lmb_edge_id as usize) {
                let mut q = m.get(&vk_symbol!("q_")).unwrap().to_owned();
                // The defining equation is k_new(i)=q_old(LMB edge i).
                // Retain this direction before the existing solve computes
                // the opposite direction used to rewrite the propagators.
                parent_coordinates.push(q.clone());
                for (src, (trgt, (_trgt_symbol, _trgt_rotated_symbol))) in
                    mom_vecs_to_symbols.iter()
                {
                    q = q.replace(src).with(trgt.clone());
                }
                q -= Atom::var(mom_vecs_to_symbols[i_lmb].1.1.1);
                system.push(q);
            }
        }
        // println!(
        //     "system: {:?}",
        //     system.iter().map(|a| a.to_string()).collect::<Vec<_>>()
        // );
        let variables = mom_vecs_to_symbols
            .iter()
            .map(|(_src, (_trgt, (trgt_symbol, _trgt_rotated_symbol)))| Atom::var(*trgt_symbol))
            .collect::<Vec<_>>();
        // println!(
        //     "variables: {:?}",
        //     variables.iter().map(|a| a.to_string()).collect::<Vec<_>>()
        // );
        let basis_change = Arc::new(Self::solve_loop_coordinates(&system, &variables)?);
        // println!(
        //     "basis_change: {:?}",
        //     basis_change
        //         .iter()
        //         .map(|a| a.to_string())
        //         .collect::<Vec<_>>()
        // );

        let mut rotated_result = input_canonical_expression.to_owned();
        for prop_id in prop_ids {
            rotated_result = rotated_result
                .replace(
                    vk_parse!(format!("prop({},edges_,q_,mUVsq_,pow_)", prop_id).as_str())
                        .unwrap()
                        .to_pattern(),
                )
                .with_map({
                    let mom_vecs_to_symbols = mom_vecs_to_symbols.clone();
                    let basis_change = basis_change.clone();
                    move |match_in| {
                        let mut q = match_in.get(vk_symbol!("q_")).unwrap().to_atom();
                        for ((src, (_trgt, _trgt_symbol)), rotated_expr) in
                            mom_vecs_to_symbols.iter().zip(basis_change.iter())
                        {
                            q = q.replace(src).with(rotated_expr.to_pattern());
                        }
                        // Map back symbols to k(i) atoms
                        for (src, (_trgt, (_trgt_symbol, trgt_rotated_symbol))) in
                            mom_vecs_to_symbols.iter()
                        {
                            q = q
                                .replace(Atom::var(*trgt_rotated_symbol).to_pattern())
                                .with(src.clone());
                        }
                        q = q.expand();
                        function!(
                            vk_symbol!("prop"),
                            Atom::num(prop_id as i64),
                            match_in.get(vk_symbol!("edges_")).unwrap().to_atom(),
                            q,
                            match_in.get(vk_symbol!("mUVsq_")).unwrap().to_atom(),
                            match_in.get(vk_symbol!("pow_")).unwrap().to_atom()
                        )
                    }
                });
        }

        // println!("Rotated topology: {}", rotated_result);

        Ok((rotated_result, parent_coordinates.into()))
    }

    pub fn get_integral(&self) -> &Integral {
        match self {
            Topology::OneLoop(i)
            | Topology::TwoLoop(i)
            | Topology::ThreeLoop(i)
            | Topology::FourLoop(i)
            | Topology::FiveLoop(i)
            | Topology::Unknown(i) => i,
        }
    }

    pub fn get_integral_mut(&mut self) -> &mut Integral {
        match self {
            Topology::OneLoop(i)
            | Topology::TwoLoop(i)
            | Topology::ThreeLoop(i)
            | Topology::FourLoop(i)
            | Topology::FiveLoop(i)
            | Topology::Unknown(i) => i,
        }
    }

    pub fn to_canonical(
        &self,
        integral: AtomView,
        replacement_rules: &ReplacementRules,
        short_form: bool,
    ) -> Atom {
        match self {
            Topology::Unknown(_) => integral
                .replace(vk_parse!("topo(props_)").unwrap().to_pattern())
                .with(vk_parse!("topo(UNKNOWN(props_))").unwrap().to_pattern()),
            t => t.get_integral().to_canonical(replacement_rules, short_form),
        }
    }

    fn match_topology_to_user_input(
        &self,
        input: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        match self {
            Topology::Unknown(unknown_integral) => {
                // println!(
                //     "\n>>>Trying to match with the unknown integral pattern {}\n<<<",
                //     self
                // );
                let undirected_input = input
                    .replace(vk_parse!("edge(x_,y_)").unwrap().to_pattern())
                    .with(vk_parse!("uedge(x_,y_)").unwrap().to_pattern());
                if undirected_input
                    .pattern_match(
                        &unknown_integral.generic_pattern.pattern,
                        Some(&unknown_integral.generic_pattern.conditions),
                        Some(&unknown_integral.generic_pattern.match_settings),
                    )
                    .next()
                    .is_some()
                {
                    // Make sure all propagators have IDs starting at 1 until n_props
                    let mut test_input = input.to_owned();
                    let mut prop_id = 1;
                    let mut prop_pattern = vk_parse!(format!("prop({},args__)", prop_id).as_str())
                        .unwrap()
                        .to_pattern();
                    while test_input
                        .pattern_match(&prop_pattern, Some(&Condition::default()), None)
                        .next()
                        .is_some()
                    {
                        test_input = test_input
                            .replace(&prop_pattern)
                            .with(Atom::num(1).to_pattern());
                        prop_id += 1;
                        prop_pattern = vk_parse!(format!("prop({},args__)", prop_id).as_str())
                            .unwrap()
                            .to_pattern();
                    }
                    test_input = test_input
                        .replace(vk_parse!("UNKNOWN(args__)").unwrap().to_pattern())
                        .with(vk_parse!("args__").unwrap().to_pattern());
                    if test_input != vk_parse!("topo(1)").unwrap() {
                        return Err(VakintError::InvalidIntegralFormat(format!(
                            "UNKNOWN integrals must have propagator ids ranging from 1 to their maximal number of propagators: {}",
                            test_input
                        )));
                    }

                    let mut loop_momenta_ids = HashSet::new();
                    for i_prop in 1..=prop_id - 1 {
                        for (mom_symbol, (_, mom_id)) in get_individual_momenta(
                            get_prop_with_id(input, i_prop)
                                .unwrap()
                                .get(&vk_symbol!("q_"))
                                .unwrap()
                                .as_view(),
                        )?
                        .iter()
                        {
                            if *mom_symbol != S.k {
                                return Err(VakintError::InvalidIntegralFormat(format!(
                                    "Unknown integrals must have only loop momenta with symbol {}: {}",
                                    S.k, mom_symbol
                                )));
                            }
                            if *mom_id < 1 {
                                return Err(VakintError::InvalidIntegralFormat(format!(
                                    "Unknown integrals must have loop momenta ids starting from 1, not {}",
                                    mom_id
                                )));
                            }
                            loop_momenta_ids.insert(*mom_id);
                        }
                    }
                    let mut sorted_loop_momenta_ids = loop_momenta_ids
                        .iter()
                        .map(|mom_id| *mom_id as usize)
                        .collect::<Vec<_>>();
                    sorted_loop_momenta_ids.sort();
                    if sorted_loop_momenta_ids
                        != (1..=sorted_loop_momenta_ids.len()).collect::<Vec<_>>()
                    {
                        return Err(VakintError::InvalidIntegralFormat(format!(
                            "Unknown integrals must have loop momenta ids ranging from 1 to their maximal number of loop momenta, not {:?}",
                            sorted_loop_momenta_ids
                        )));
                    }
                    let n_loops = *sorted_loop_momenta_ids.last().unwrap();
                    //println!("Found a match!");
                    let mut replacement_rules = ReplacementRules {
                        canonical_topology: Topology::Unknown(unknown_integral.clone()),
                        ..ReplacementRules::default()
                    };
                    replacement_rules.canonical_expression_substitutions.insert(
                        vk_parse!("integral").unwrap(),
                        input
                            .replace(vk_parse!("UNKNOWN(int_)").unwrap().to_pattern())
                            .with(vk_parse!("int_").unwrap().to_pattern()),
                    );
                    replacement_rules.canonical_expression_substitutions.insert(
                        vk_parse!("n_props").unwrap(),
                        Atom::num((prop_id - 1) as i64),
                    );
                    replacement_rules
                        .canonical_expression_substitutions
                        .insert(vk_parse!("n_loops").unwrap(), Atom::num(n_loops as i64));
                    Ok(Some(replacement_rules))
                } else {
                    //println!("Does not match!");
                    Ok(None)
                }
            }
            t => {
                // println!("\n>>>Trying to match: {}\n<<<", self);
                if let Some(mut replacement_rule) =
                    t.get_integral().match_integral_to_user_input(input)?
                {
                    replacement_rule.canonical_topology = t.clone();
                    // println!("Found a match! ->\n{}", replacement_rule);
                    Ok(Some(replacement_rule))
                } else {
                    Ok(None)
                }
            }
        }
    }
}
