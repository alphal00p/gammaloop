pub mod diagram_generator;

use ahash::{AHashMap, HashMap};
use smartstring::{LazyCompact, SmartString};
use std::{fmt, str::FromStr};
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
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum NumeratorAwareGraphGroupingOption {
    NoGrouping,
    OnlyDetectZeroes,
    GroupIdenticalGraphUpToSign,
    GroupIdenticalGraphUpToScalarRescaling,
}

impl fmt::Display for NumeratorAwareGraphGroupingOption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::NoGrouping => "No grouping",
                Self::OnlyDetectZeroes => "Only detect zero numerators",
                Self::GroupIdenticalGraphUpToSign => "Group identical graphs up to a sign",
                Self::GroupIdenticalGraphUpToScalarRescaling => {
                    "Group identical graphs up to a scalar rescaling"
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
            "group_identical_graphs_up_to_sign" => Ok(Self::GroupIdenticalGraphUpToSign),
            "group_identical_graphs_up_to_scalar_rescaling" => {
                Ok(Self::GroupIdenticalGraphUpToScalarRescaling)
            }
            _ => Err(FeynGenError::GenericError(format!(
                "Invalid grouping option: {}",
                s
            ))),
        }
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
    graph: &SymbolicaGraph<NodeColor, &str>,
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

    pub fn allow_tadpoles(&self) -> bool {
        !self
            .0
            .iter()
            .any(|f| matches!(f, FeynGenFilter::TadpolesFilter(_)))
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

    pub fn get_coupling_orders(&self) -> Option<&HashMap<String, usize>> {
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

    #[allow(clippy::type_complexity)]
    pub fn apply_filters<NodeColor: diagram_generator::NodeColorFunctions>(
        &self,
        graphs: &mut Vec<(SymbolicaGraph<NodeColor, &str>, String)>,
    ) -> Result<(), FeynGenError> {
        for filter in self.0.iter() {
            #[allow(clippy::single_match)]
            match filter {
                FeynGenFilter::CouplingOrders(orders) => {
                    graphs.retain(|(g, _)| {
                        let graph_coupling_orders = get_coupling_orders(g);
                        orders.iter().all(|(k, v)| {
                            graph_coupling_orders
                                .get(&SmartString::from(k))
                                .map_or(0 == *v, |o| *o == *v)
                        })
                    });
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

#[derive(Debug, Clone)]
pub enum FeynGenFilter {
    SelfEnergyFilter(SelfEnergyFilterOptions),
    TadpolesFilter(TadpolesFilterOptions),
    ZeroSnailsFilter(SnailFilterOptions),
    ParticleVeto(Vec<i64>),
    MaxNumberOfBridges(usize),
    CouplingOrders(HashMap<String, usize>),
    LoopCountRange((usize, usize)),
    PerturbativeOrders(HashMap<String, usize>),
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
                Self::MaxNumberOfBridges(n) => format!("MaxNumberOfBridges({})", n),
                Self::TadpolesFilter(opts) => format!("NoTadpoles({})", opts),
                Self::ZeroSnailsFilter(opts) => format!("NoZeroSnails({})", opts),
                Self::CouplingOrders(orders) => format!(
                    "CouplingOrders({})",
                    orders
                        .iter()
                        .map(|(k, v)| format!("{}={}", k, v))
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
            }
        )
    }
}

#[derive(Debug, Clone)]
pub struct FeynGenOptions {
    pub generation_type: GenerationType,
    pub initial_pdgs: Vec<i64>,
    pub final_pdgs: Vec<i64>,
    pub loop_count_range: (usize, usize),
    pub symmetrize_initial_states: bool,
    pub symmetrize_final_states: bool,
    pub symmetrize_left_right_states: bool,
    pub amplitude_filters: FeynGenFilters,
    pub cross_section_filters: FeynGenFilters,
}

impl fmt::Display for FeynGenOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Generation type: {}{}\nInitial PDGs: {:?}{}\nFinal PDGs: {:?}{}\nLoop count: {}\nAmplitude filters:{}{}\nCross-section filters:{}{}",
            self.generation_type,
            if self.symmetrize_left_right_states { " (left-right symmetrized)" } else { "" },
            self.initial_pdgs,
            if self.symmetrize_initial_states { " (symmetrized)" } else { "" },
            self.final_pdgs,
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
