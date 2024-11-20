pub mod diagram_generator;

use ahash::HashMap;
use smartstring::{LazyCompact, SmartString};
use std::{fmt, str::FromStr};
use symbolica::graph::Graph as SymbolicaGraph;
use thiserror::Error;

use crate::model::Model;

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

pub fn get_coupling_orders(
    model: &Model,
    graph: &SymbolicaGraph<(i32, SmartString<LazyCompact>), &str>,
) -> HashMap<String, usize> {
    let mut coupling_orders = HashMap::default();
    for node in graph.nodes() {
        if node.data.1 == "external" {
            continue;
        }
        let mut node_coupling_orders: HashMap<String, usize> = HashMap::default();
        model
            .get_vertex_rule(&node.data.1)
            .couplings
            .iter()
            .for_each(|cs| {
                cs.iter().for_each(|c_opt| {
                    if let Some(c) = c_opt {
                        c.orders.iter().for_each(|(coupling_order, &weight)| {
                            if let Some(curr_order) =
                                node_coupling_orders.get_mut(&coupling_order.clone().to_string())
                            {
                                if *curr_order < weight {
                                    *curr_order = weight;
                                }
                            } else {
                                node_coupling_orders.insert(coupling_order.clone().into(), weight);
                            }
                        });
                    }
                })
            });
        for (k, v) in node_coupling_orders {
            *coupling_orders.entry(k).or_insert_with(|| 0) += v;
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

    pub fn get_particle_vetos(&self) -> Option<&FeynGenFilter> {
        self.0
            .iter()
            .find(|f| matches!(f, FeynGenFilter::ParticleVeto(_)))
    }

    #[allow(clippy::type_complexity)]
    pub fn apply_filters(
        &self,
        model: &Model,
        graphs: &mut Vec<(
            SymbolicaGraph<(i32, SmartString<LazyCompact>), &str>,
            String,
        )>,
    ) -> Result<(), FeynGenError> {
        for filter in self.0.iter() {
            #[allow(clippy::single_match)]
            match filter {
                FeynGenFilter::CouplingOrders(orders) => {
                    graphs.retain(|(g, _)| {
                        let graph_coupling_orders = get_coupling_orders(model, g);
                        orders.iter().all(|(k, v)| {
                            graph_coupling_orders.get(k).map_or(0 == *v, |o| *o == *v)
                        })
                    });
                }
                FeynGenFilter::MaxNumberOfBridges(_)
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
    pub filters: FeynGenFilters,
}

impl fmt::Display for FeynGenOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Generation type: {}{}\nInitial PDGs: {:?}{}\nFinal PDGs: {:?}{}\nLoop count: {}\nFilters:{}{}",
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
            if self.filters.0.is_empty() { " None" } else {"\n"},
            if self.filters.0.is_empty() { "".into() } else { self.filters
                .0
                .iter()
                .map(|f| format!(" > {}", f))
                .collect::<Vec<String>>()
                .join("\n")
            }
        )
    }
}
