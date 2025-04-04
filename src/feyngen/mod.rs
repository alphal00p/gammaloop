pub mod diagram_generator;

use ahash::{AHashMap, HashMap};
use diagram_generator::EdgeColor;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use itertools::Itertools;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use smartstring::{LazyCompact, SmartString};
use std::ops::RangeInclusive;
use std::{fmt, str::FromStr};
use symbolica::atom::Atom;
use symbolica::graph::Graph as SymbolicaGraph;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FeynGenError {
    #[error("{0}")]
    GenericError(String),
    #[error("Invalid loop momentum basis | {0}")]
    LoopMomentumBasisError(String),
    #[error("Could not convert symbolica graph symmetry factor to an integer: {0}")]
    SymmetryFactorError(String),
    #[error("Could not numerically evaluate numerator: {0}")]
    NumeratorEvaluationError(String),
}

#[derive(Debug, Clone, Eq, PartialEq)]

pub struct GraphGroupingOptions {
    pub numerical_sample_seed: u16,
    pub number_of_numerical_samples: usize,
    pub differentiate_particle_masses_only: bool,
    pub fully_numerical_substitution_when_comparing_numerators: bool,
    pub test_canonized_numerator: bool,
}

impl Default for GraphGroupingOptions {
    fn default() -> Self {
        Self {
            numerical_sample_seed: 3,
            number_of_numerical_samples: 5,
            differentiate_particle_masses_only: true,
            fully_numerical_substitution_when_comparing_numerators: true,
            test_canonized_numerator: false,
        }
    }
}

impl fmt::Display for GraphGroupingOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "differentiate_masses_only={}, test_canonized_numerator={}, #samples={}, seed={}",
            self.numerical_sample_seed,
            self.number_of_numerical_samples,
            self.differentiate_particle_masses_only,
            self.test_canonized_numerator
        )
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum NumeratorAwareGraphGroupingOption {
    NoGrouping,
    OnlyDetectZeroes,
    GroupIdenticalGraphUpToSign(GraphGroupingOptions),
    GroupIdenticalGraphUpToScalarRescaling(GraphGroupingOptions),
}

impl fmt::Display for NumeratorAwareGraphGroupingOption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::NoGrouping => "no grouping",
                Self::OnlyDetectZeroes => "only detect zero numerators",
                Self::GroupIdenticalGraphUpToSign(_opts) => "up to a sign",
                Self::GroupIdenticalGraphUpToScalarRescaling(_opts) => {
                    "up to a scalar rescaling"
                }
            }
        )
    }
}

impl FromStr for NumeratorAwareGraphGroupingOption {
    type Err = FeynGenError;

    fn from_str(s: &str) -> Result<Self, FeynGenError> {
        match s {
            "no_grouping" => Ok(Self::NoGrouping),
            "only_detect_zeroes" => Ok(Self::OnlyDetectZeroes),
            "group_identical_graphs_up_to_sign" => Ok(Self::GroupIdenticalGraphUpToSign(
                GraphGroupingOptions::default(),
            )),
            "group_identical_graphs_up_to_scalar_rescaling" => Ok(
                Self::GroupIdenticalGraphUpToScalarRescaling(GraphGroupingOptions::default()),
            ),
            _ => Err(FeynGenError::GenericError(format!(
                "Invalid grouping option: {}",
                s
            ))),
        }
    }
}

impl NumeratorAwareGraphGroupingOption {
    pub fn set_options(&mut self) -> Option<&mut GraphGroupingOptions> {
        match self {
            Self::NoGrouping => None,
            Self::OnlyDetectZeroes => None,
            Self::GroupIdenticalGraphUpToSign(opts) => Some(opts),
            Self::GroupIdenticalGraphUpToScalarRescaling(opts) => Some(opts),
        }
    }

    pub fn get_options(&self) -> Option<&GraphGroupingOptions> {
        match self {
            Self::NoGrouping => None,
            Self::OnlyDetectZeroes => None,
            Self::GroupIdenticalGraphUpToSign(opts) => Some(opts),
            Self::GroupIdenticalGraphUpToScalarRescaling(opts) => Some(opts),
        }
    }

    pub fn description(&self) -> String {
        format!(
            "{}{}",
            self,
            self.get_options().map_or("".into(), |o| format!("({})", o))
        )
    }

    pub fn new_with_attributes(
        strategy: &str,
        seed: Option<u16>,
        num_samples: Option<usize>,
        differentiate_particle_masses_only: Option<bool>,
        fully_numerical_substitution_when_comparing_numerators: Option<bool>,
        test_canonized_numerator: Option<bool>,
    ) -> Result<Self, FeynGenError> {
        let mut opt = NumeratorAwareGraphGroupingOption::from_str(strategy)?;
        let grouping_options = opt.set_options();
        if let Some(grouping_options) = grouping_options {
            if let Some(seed) = seed {
                grouping_options.numerical_sample_seed = seed;
            }
            if let Some(num_samples) = num_samples {
                grouping_options.number_of_numerical_samples = num_samples;
            }
            if let Some(differentiate_particle_masses_only) = differentiate_particle_masses_only {
                grouping_options.differentiate_particle_masses_only =
                    differentiate_particle_masses_only;
            }
            if let Some(test_canonized_numerator) = test_canonized_numerator {
                grouping_options.test_canonized_numerator = test_canonized_numerator;
            }
            if let Some(fully_numerical_substitution_when_comparing_numerators) =
                fully_numerical_substitution_when_comparing_numerators
            {
                grouping_options.fully_numerical_substitution_when_comparing_numerators =
                    fully_numerical_substitution_when_comparing_numerators;
            }
        }

        Ok(opt)
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum GenerationType {
    Amplitude,
    CrossSection,
}

impl fmt::Display for GenerationType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Amplitude => "Amplitude",
                Self::CrossSection => "Cross-section",
            }
        )
    }
}

impl FromStr for GenerationType {
    type Err = FeynGenError;

    fn from_str(s: &str) -> Result<Self, FeynGenError> {
        match s {
            "amplitude" => Ok(Self::Amplitude),
            "cross_section" => Ok(Self::CrossSection),
            _ => Err(FeynGenError::GenericError(format!(
                "Invalid generation type: {}",
                s
            ))),
        }
    }
}

pub fn get_coupling_orders<NodeColor: diagram_generator::NodeColorFunctions>(
    graph: &SymbolicaGraph<NodeColor, EdgeColor>,
) -> AHashMap<SmartString<LazyCompact>, usize> {
    let mut coupling_orders = AHashMap::default();
    for node in graph.nodes() {
        for (k, v) in node.data.coupling_orders() {
            *coupling_orders.entry(k).or_insert(0) += v;
        }
    }
    coupling_orders
}

#[derive(Debug, Clone)]
pub struct FeynGenFilters(pub Vec<FeynGenFilter>);

impl FeynGenFilters {
    pub fn get_max_bridge(&self) -> Option<usize> {
        self.0.iter().find_map(|f| match f {
            FeynGenFilter::MaxNumberOfBridges(n) => Some(*n),
            _ => None,
        })
    }

    pub fn get_blob_range(&self) -> Option<&RangeInclusive<usize>> {
        self.0.iter().find_map(|f| {
            if let FeynGenFilter::BlobRange(v) = f {
                Some(v)
            } else {
                None
            }
        })
    }

    pub fn get_spectator_range(&self) -> Option<&RangeInclusive<usize>> {
        self.0.iter().find_map(|f| {
            if let FeynGenFilter::SpectatorRange(v) = f {
                Some(v)
            } else {
                None
            }
        })
    }

    pub fn allow_tadpoles(&self) -> bool {
        !self
            .0
            .iter()
            .any(|f| matches!(f, FeynGenFilter::TadpolesFilter(_)))
    }

    pub fn filter_cross_section_tadpoles(&self) -> bool {
        self.0.iter().any(|f| {
            matches!(
                f,
                FeynGenFilter::SewedFilter(SewedFilterOptions {
                    filter_tadpoles: true,
                    ..
                })
            )
        })
    }

    pub fn get_particle_vetos(&self) -> Option<&[i64]> {
        self.0.iter().find_map(|f| {
            if let FeynGenFilter::ParticleVeto(v) = f {
                Some(v.as_slice())
            } else {
                None
            }
        })
    }

    pub fn get_coupling_orders(&self) -> Option<&HashMap<String, (usize, Option<usize>)>> {
        self.0.iter().find_map(|f| {
            if let FeynGenFilter::CouplingOrders(o) = f {
                Some(o)
            } else {
                None
            }
        })
    }

    pub fn get_perturbative_orders(&self) -> Option<&HashMap<String, usize>> {
        self.0.iter().find_map(|f| {
            if let FeynGenFilter::PerturbativeOrders(o) = f {
                Some(o)
            } else {
                None
            }
        })
    }

    pub fn get_loop_count_range(&self) -> Option<(usize, usize)> {
        self.0.iter().find_map(|f| {
            if let FeynGenFilter::LoopCountRange(o) = f {
                Some(*o)
            } else {
                None
            }
        })
    }

    pub fn get_fermion_loop_count_range(&self) -> Option<(usize, usize)> {
        self.0.iter().find_map(|f: &FeynGenFilter| {
            if let FeynGenFilter::FermionLoopCountRange(o) = f {
                Some(*o)
            } else {
                None
            }
        })
    }

    #[allow(clippy::type_complexity)]
    pub fn apply_filters<NodeColor: diagram_generator::NodeColorFunctions + Send + Sync + Clone>(
        &self,
        graphs: &mut Vec<(SymbolicaGraph<NodeColor, EdgeColor>, Atom)>,
        pool: &rayon::ThreadPool,
        progress_bar_style: &ProgressStyle,
    ) -> Result<(), FeynGenError> {
        for filter in self.0.iter() {
            match filter {
                FeynGenFilter::CouplingOrders(orders) => {
                    let bar = ProgressBar::new(graphs.len() as u64);
                    bar.set_style(progress_bar_style.clone());
                    bar.set_message("Applying coupling orders constraints...");
                    pool.install(|| {
                        *graphs = graphs
                            .par_iter_mut()
                            .progress_with(bar.clone())
                            .filter(|(g, _)| {
                                let graph_coupling_orders = get_coupling_orders(g);
                                let a = orders.iter().all(|(k, (v_min, v_max))| {
                                    graph_coupling_orders
                                        .get(&SmartString::from(k))
                                        .map_or(0 == *v_min, |o| {
                                            *o >= *v_min && (*v_max).is_none_or(|max| *o <= max)
                                        })
                                });
                                // if a {
                                //     info!(
                                //         "Coupling orders constraints satisfied for graph {}",
                                //         g.to_dot()
                                //     );
                                //     info!("{:?}", graph_coupling_orders);
                                // }
                                a
                            })
                            .map(|(g, sf)| (g.clone(), sf.clone()))
                            .collect::<Vec<_>>()
                    });
                    bar.finish_and_clear();
                }
                FeynGenFilter::LoopCountRange((loop_count_min, loop_count_max)) => {
                    graphs.retain(|(g, _)| {
                        g.num_loops() >= *loop_count_min && g.num_loops() <= *loop_count_max
                    });
                }
                FeynGenFilter::PerturbativeOrders(_)
                | FeynGenFilter::MaxNumberOfBridges(_)
                | FeynGenFilter::SelfEnergyFilter(_)
                | FeynGenFilter::TadpolesFilter(_)
                | FeynGenFilter::ZeroSnailsFilter(_)
                | FeynGenFilter::FermionLoopCountRange(_)
                | FeynGenFilter::SewedFilter(_)
                | FeynGenFilter::FactorizedLoopTopologiesCountRange(_)
                | FeynGenFilter::BlobRange(_)
                | FeynGenFilter::SpectatorRange(_)
                | FeynGenFilter::ParticleVeto(_) => {} // These other filters are implemented directly during diagram generation
            }
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SelfEnergyFilterOptions {
    pub veto_self_energy_of_massive_lines: bool,
    pub veto_self_energy_of_massless_lines: bool,
    pub veto_only_scaleless_self_energy: bool,
}

impl Default for SelfEnergyFilterOptions {
    fn default() -> Self {
        Self {
            veto_self_energy_of_massive_lines: true,
            veto_self_energy_of_massless_lines: true,
            veto_only_scaleless_self_energy: false,
        }
    }
}

impl fmt::Display for SelfEnergyFilterOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut descr = vec![];
        if self.veto_self_energy_of_massive_lines && !self.veto_self_energy_of_massless_lines {
            descr.push("only massive legs")
        } else if !self.veto_self_energy_of_massive_lines && self.veto_self_energy_of_massless_lines
        {
            descr.push("only massless legs")
        };
        if self.veto_only_scaleless_self_energy {
            descr.push("only scaleless self-energies")
        };
        write!(f, "{}", descr.join(" | "))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SnailFilterOptions {
    pub veto_snails_attached_to_massive_lines: bool,
    pub veto_snails_attached_to_massless_lines: bool,
    pub veto_only_scaleless_snails: bool,
}

impl Default for SnailFilterOptions {
    fn default() -> Self {
        Self {
            veto_snails_attached_to_massive_lines: false,
            veto_snails_attached_to_massless_lines: true,
            veto_only_scaleless_snails: false,
        }
    }
}

impl fmt::Display for SnailFilterOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut descr = vec![];
        if self.veto_snails_attached_to_massive_lines
            && !self.veto_snails_attached_to_massless_lines
        {
            descr.push("only attached to massive legs")
        } else if !self.veto_snails_attached_to_massive_lines
            && self.veto_snails_attached_to_massless_lines
        {
            descr.push("only attached to massless legs")
        };
        if self.veto_only_scaleless_snails {
            descr.push("only scaleless snails")
        };
        write!(f, "{}", descr.join(" | "))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct TadpolesFilterOptions {
    pub veto_tadpoles_attached_to_massive_lines: bool,
    pub veto_tadpoles_attached_to_massless_lines: bool,
    pub veto_only_scaleless_tadpoles: bool,
}

impl Default for TadpolesFilterOptions {
    fn default() -> Self {
        Self {
            veto_tadpoles_attached_to_massive_lines: true,
            veto_tadpoles_attached_to_massless_lines: true,
            veto_only_scaleless_tadpoles: false,
        }
    }
}

impl fmt::Display for TadpolesFilterOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut descr = vec![];
        if self.veto_tadpoles_attached_to_massive_lines
            && !self.veto_tadpoles_attached_to_massless_lines
        {
            descr.push("only attached to massive legs")
        } else if !self.veto_tadpoles_attached_to_massive_lines
            && self.veto_tadpoles_attached_to_massless_lines
        {
            descr.push("only attached to massless legs")
        };
        if self.veto_only_scaleless_tadpoles {
            descr.push("only scaleless tadpoles")
        };
        write!(f, "{}", descr.join(" | "))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SewedFilterOptions {
    pub filter_tadpoles: bool,
}

#[derive(Debug, Clone)]
pub enum FeynGenFilter {
    SelfEnergyFilter(SelfEnergyFilterOptions),
    TadpolesFilter(TadpolesFilterOptions),
    ZeroSnailsFilter(SnailFilterOptions),
    SewedFilter(SewedFilterOptions),
    /// A list of vetoed pdgs
    ParticleVeto(Vec<i64>),
    MaxNumberOfBridges(usize),
    /// A map between the coupling order name and a range of orders, inclusive, with an optional upper bound
    CouplingOrders(HashMap<String, (usize, Option<usize>)>),
    /// A range of loop counts, inclusive
    LoopCountRange((usize, usize)),
    /// A range of blob counts, inclusive
    BlobRange(RangeInclusive<usize>),
    SpectatorRange(RangeInclusive<usize>),
    PerturbativeOrders(HashMap<String, usize>),
    FermionLoopCountRange((usize, usize)),
    FactorizedLoopTopologiesCountRange((usize, usize)),
}

impl fmt::Display for FeynGenFilter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::SelfEnergyFilter(opts) => format!("NoExternalSelfEnergy({})", opts),
                Self::ParticleVeto(pdgs) => format!(
                    "ParticleVeto({})",
                    pdgs.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join("|")
                ),
                Self::SpectatorRange(r) => format!("SpectatorRange({:?})", r),
                Self::BlobRange(r) => format!("BlobRange({:?})", r),
                Self::MaxNumberOfBridges(n) => format!("MaxNumberOfBridges({})", n),
                Self::TadpolesFilter(opts) => format!("NoTadpoles({})", opts),
                Self::ZeroSnailsFilter(opts) => format!("NoZeroSnails({})", opts),
                Self::CouplingOrders(orders) => format!(
                    "CouplingOrders({})",
                    orders
                        .iter()
                        .map(|(k, (v_min, v_max_opt))| {
                            if let Some(v_max) = v_max_opt {
                                if v_min == v_max {
                                    format!("{}=={}", k, v_min)
                                } else {
                                    format!("{}=[{}..{}]", k, v_min, v_max)
                                }
                            } else {
                                format!("{}>={}", k, v_min)
                            }
                        })
                        .collect::<Vec<String>>()
                        .join("|")
                ),
                Self::PerturbativeOrders(orders) => format!(
                    "PerturbativeOrders({})",
                    orders
                        .iter()
                        .map(|(k, v)| format!("{}={}", k, v))
                        .collect::<Vec<String>>()
                        .join("|")
                ),
                Self::LoopCountRange((loop_count_min, loop_count_max)) =>
                    format!("LoopCountRange({{{},{}}})", loop_count_min, loop_count_max),
                Self::FermionLoopCountRange((loop_count_min, loop_count_max)) => format!(
                    "FermionLoopCountRange({{{},{}}})",
                    loop_count_min, loop_count_max
                ),
                Self::FactorizedLoopTopologiesCountRange((loop_count_min, loop_count_max)) =>
                    format!(
                        "NFactorizableLoopRange({{{},{}}})",
                        loop_count_min, loop_count_max
                    ),
                Self::SewedFilter(SewedFilterOptions { filter_tadpoles }) => format!(
                    "SewedCrossSectionFilter(filter_tadpoles={{{}}})",
                    filter_tadpoles
                ),
            }
        )
    }
}

#[derive(Debug, Clone)]
pub struct FeynGenOptions {
    pub generation_type: GenerationType,
    pub initial_pdgs: Vec<i64>,
    pub final_pdgs_lists: Vec<Vec<i64>>,
    pub loop_count_range: (usize, usize),
    pub symmetrize_initial_states: bool,
    pub symmetrize_final_states: bool,
    pub symmetrize_left_right_states: bool,
    pub allow_symmetrization_of_external_fermions_in_amplitudes: bool,
    pub max_multiplicity_for_fast_cut_filter: usize,
    pub amplitude_filters: FeynGenFilters,
    pub cross_section_filters: FeynGenFilters,
}

impl fmt::Display for FeynGenOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Generation type: {}{}{}\nInitial PDGs: {:?}{}\nFinal PDGs: {}{}\nLoop count: {}\nAmplitude filters:{}{}\nCross-section filters:{}{}",
            self.generation_type,
            if self.symmetrize_left_right_states { " (left-right symmetrized)" } else { "" },
            if self.allow_symmetrization_of_external_fermions_in_amplitudes
                && self.generation_type == GenerationType::Amplitude
                && (self.symmetrize_initial_states  || self.symmetrize_final_states || self.symmetrize_left_right_states)
                { " (allowing fermion symmetrization)" } else { "" },
            self.initial_pdgs,
            if self.symmetrize_initial_states { " (symmetrized)" } else { "" },
            if self.final_pdgs_lists.len() == 1 {
                format!("{:?}",self.final_pdgs_lists[0])
            } else {
                format!("[ {} ]", self.final_pdgs_lists.iter().map(|pdgs| format!("{:?}", pdgs)).join(" | "))
            },
            if self.symmetrize_final_states { " (symmetrized)" } else { "" },
            if self.loop_count_range.0 == self.loop_count_range.1 {
                format!("{}", self.loop_count_range.0)
            } else {
                format!("{:?}", self.loop_count_range)
            },
            if self.amplitude_filters.0.is_empty() { " None" } else {"\n"},
            if self.amplitude_filters.0.is_empty() { "".into() } else { self.amplitude_filters
                .0
                .iter()
                .map(|f| format!(" > {}", f))
                .collect::<Vec<String>>()
                .join("\n")
            },
            if self.cross_section_filters.0.is_empty() { " None" } else {"\n"},
            if self.cross_section_filters.0.is_empty() { "".into() } else { self.cross_section_filters
                .0
                .iter()
                .map(|f| format!(" > {}", f))
                .collect::<Vec<String>>()
                .join("\n")
            }
        )
    }
}

pub mod half_edge_filters;
#[cfg(test)]
pub mod test;
