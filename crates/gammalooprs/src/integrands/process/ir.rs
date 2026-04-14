use std::{collections::HashSet, fmt::Display};

use color_eyre::eyre::Result;
use colored::Colorize;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};
use nalgebra::DVector;
use rand::Rng;
use spenso::algebra::complex::Complex;
use symbolica::{
    domains::float::{Real, RealLike},
    numerical_integration::{MonteCarloRng, Sample},
};
use tabled::{builder::Builder, settings::Style};
use tracing::warn;

use crate::{
    DependentMomentaConstructor,
    graph::{FeynmanGraph, LmbError, lmb::LMBwithEdges},
    integrands::{
        HasIntegrand,
        evaluation::EvaluationResult,
        process::{
            ProcessIntegrandImpl,
            amplitude::{AmplitudeGraphTerm, AmplitudeIntegrand},
            cross_section::{CrossSectionGraphTerm, CrossSectionIntegrand},
        },
    },
    model::Model,
    momentum::{
        ThreeMomentum,
        sample::{LoopIndex, LoopMomenta, MomentumSample},
    },
    settings::{
        RuntimeSettings, SamplingSettings,
        runtime::{
            DiscreteGraphSamplingSettings, DiscreteGraphSamplingType, ParameterizationMapping,
            ParameterizationMode, ParameterizationSettings,
        },
    },
    utils::{F, FloatLike, box_muller},
    uv::profile::logspace,
};
use varpro::{
    prelude::SeparableModelBuilder, problem::SeparableProblemBuilder, solvers::levmar::LevMarSolver,
};

/// The range is from 10^start to 10^end.
pub struct IRProfileSetting {
    pub lambda_exp_start: f64,
    pub lambda_exp_end: f64,
    pub steps: usize,
    pub seed: u64,
    pub select_limits_and_graphs: Option<String>,
}

impl AmplitudeGraphTerm {
    fn enumerate_ir_limits(&self) -> Vec<IrLimit> {
        let mut limits: HashSet<IrLimit> = HashSet::new();

        let massless_edges: Vec<EdgeIndex> = self
            .graph
            .iter_edges_of(&!self.graph.tree_edges.clone())
            .filter_map(|(_a, b, c)| {
                if c.data.particle.is_massless() {
                    Some(b)
                } else {
                    None
                }
            })
            .collect();

        for subset in massless_edges.iter().powerset() {
            if subset.is_empty() {
                continue;
            }
            let _lmb = match self.lmb_with_loop_edges(subset.as_slice()) {
                Ok(lmb) => lmb,
                Err(err) => match err {
                    LmbError::NotLoopEdges { .. } => {
                        // warn!("{loop_edges} is not a valid loop edge subset");
                        continue;
                    }
                    a => panic!("Failed to build IR loop momentum basis for cut graph:\n{a}"),
                },
            };

            let ir_limit = IrLimit::new_pure_soft(subset.into_iter().copied().collect());

            limits.insert(ir_limit);
        }

        limits.into_iter().sorted().collect()
    }
}

impl CrossSectionGraphTerm {
    fn enumerate_ir_limits(&self) -> Vec<IrLimit> {
        let mut limits: HashSet<IrLimit> = HashSet::new();
        let loop_count = self.graph.loop_momentum_basis.loop_edges.len();

        for cut_group in self.raised_data.raised_cut_groups.iter() {
            let mut limits_of_cut: HashSet<IrLimit> = HashSet::new();

            let representative_cut_esurface = &self.cut_esurface[*cut_group.cuts.first().unwrap()];
            let massless_edges_in_cut = representative_cut_esurface
                .energies
                .iter()
                .filter(|edge_id| self.graph[**edge_id].particle.is_massless())
                .copied()
                .collect_vec();

            if massless_edges_in_cut.len() >= 2 {
                let subsets = massless_edges_in_cut
                    .iter()
                    .powerset()
                    .filter(|subset| {
                        subset.len() >= 2
                            && subset.len() <= loop_count
                            && subset.len() < representative_cut_esurface.energies.len()
                    })
                    .collect_vec();

                for subset in subsets {
                    let ir_limit =
                        IrLimit::new_pure_colinear(subset.into_iter().copied().collect());
                    limits_of_cut.insert(ir_limit);
                }
            }

            if !massless_edges_in_cut.is_empty() {
                let subsets = massless_edges_in_cut
                    .iter()
                    .powerset()
                    .filter(|subset| {
                        !subset.is_empty()
                            && subset.len() <= loop_count
                            && subset.len() < representative_cut_esurface.energies.len()
                    })
                    .collect_vec();
                for subset in subsets {
                    let ir_limit = IrLimit::new_pure_soft(subset.into_iter().copied().collect());
                    limits_of_cut.insert(ir_limit);
                }
            }

            for limit in limits_of_cut.drain() {
                limits.insert(limit);
            }
        }

        limits.into_iter().sorted().collect()
    }
}

pub struct IrLimitTestReport {
    pub all_passed: bool,
    pub results_per_graph: Vec<GraphIRLimitReport>,
}

pub struct GraphIRLimitReport {
    pub graph_name: String,
    pub all_limits_passed: bool,
    pub single_limit_reports: Vec<SingleLimitReport>,
}

pub struct SingleLimitReport {
    pub limit_name: String,
    pub passed: bool,
    pub power_law_fit: PowerLawFit,
    pub scaling: f64,
    num_soft: usize,
}

impl Display for IrLimitTestReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let overall_status = if self.all_passed {
            "PASS".green().bold()
        } else {
            "FAIL".red().bold()
        };

        let passed_graphs = self
            .results_per_graph
            .iter()
            .filter(|graph_report| graph_report.all_limits_passed)
            .count();

        writeln!(
            f,
            "IR limit tests: {} ({}/{})",
            overall_status,
            passed_graphs,
            self.results_per_graph.len()
        )?;

        let mut graph_summary_table = Builder::new();
        graph_summary_table.push_record(["graph", "status", "passed", "total"]);

        for graph_report in &self.results_per_graph {
            let graph_status = if graph_report.all_limits_passed {
                "PASS".green().bold().to_string()
            } else {
                "FAIL".red().bold().to_string()
            };

            let passed_limits = graph_report
                .single_limit_reports
                .iter()
                .filter(|report| report.passed)
                .count();

            graph_summary_table.push_record([
                graph_report.graph_name.clone(),
                graph_status,
                passed_limits.to_string(),
                graph_report.single_limit_reports.len().to_string(),
            ]);
        }

        writeln!(f, "{}", graph_summary_table.build().with(Style::rounded()))?;

        for graph_report in &self.results_per_graph {
            writeln!(f)?;
            writeln!(f, "{graph_report}")?;
        }

        Ok(())
    }
}

impl Display for GraphIRLimitReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let graph_status = if self.all_limits_passed {
            "PASS".green().bold()
        } else {
            "FAIL".red().bold()
        };

        let passed_limits = self
            .single_limit_reports
            .iter()
            .filter(|report| report.passed)
            .count();

        writeln!(
            f,
            "  {} {} ({}/{})",
            graph_status,
            self.graph_name.bold(),
            passed_limits,
            self.single_limit_reports.len()
        )?;

        let mut limit_table = Builder::new();
        limit_table.push_record([
            "status",
            "limit",
            "scaling",
            "p",
            "coefficient",
            "const_offset",
            "r_squared",
            "n_soft",
        ]);

        for report in &self.single_limit_reports {
            let status = if report.passed {
                "PASS".green().bold().to_string()
            } else {
                "FAIL".red().bold().to_string()
            };

            limit_table.push_record([
                status,
                report.limit_name.clone(),
                format!("{:+.4}", report.scaling),
                format!("{:+.4}", report.power_law_fit.exponent),
                format!("{:+.4e}", report.power_law_fit.coefficient),
                format!("{:+.4e}", report.power_law_fit.constant_offset),
                format!("{:.4}", report.power_law_fit.r_squared),
                report.num_soft.to_string(),
            ]);
        }

        write!(f, "{}", limit_table.build().with(Style::rounded()))?;

        Ok(())
    }
}

impl Display for SingleLimitReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let status = if self.passed {
            "PASS".green().bold()
        } else {
            "FAIL".red().bold()
        };

        write!(
            f,
            "{} {} | scaling={:+.4} | p={:+.4} | coeff={:+.4e} | const={:+.4e} | R²={:.4} | n_soft={}",
            status,
            self.limit_name,
            self.scaling,
            self.power_law_fit.exponent,
            self.power_law_fit.coefficient,
            self.power_law_fit.constant_offset,
            self.power_law_fit.r_squared,
            self.num_soft
        )
    }
}

impl AmplitudeIntegrand {
    pub fn ir_profile_completion_entries(&self) -> Vec<(String, Vec<String>)> {
        self.enumerate_ir_limits()
            .into_iter()
            .map(|(graph_name, limits)| {
                (
                    graph_name,
                    limits.into_iter().map(|limit| limit.to_string()).collect(),
                )
            })
            .collect()
    }

    pub fn test_ir(
        &mut self,
        ir_profile_settings: &IRProfileSetting,
        model: &Model,
    ) -> Result<IrLimitTestReport> {
        // override the sampling to be in momentum space
        self.settings.sampling = SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
            sample_orientations: false,
            sampling_type: DiscreteGraphSamplingType::Default(ParameterizationSettings {
                mode: ParameterizationMode::MomentumSpace,
                mapping: ParameterizationMapping::default(),
                b: 10.0,
            }),
        });

        self.warm_up(model)?;

        let mut rng = MonteCarloRng::new(ir_profile_settings.seed, 0);

        let limits_to_check =
            if let Some(select_limits_and_graphs) = &ir_profile_settings.select_limits_and_graphs {
                self.parse_select_limits_and_graphs(select_limits_and_graphs)?
            } else {
                self.enumerate_ir_limits()
            };

        let mut result = IrLimitTestReport {
            all_passed: false,
            results_per_graph: Vec::new(),
        };

        for (graph_name, limits) in limits_to_check {
            let mut graph_report = GraphIRLimitReport {
                graph_name: graph_name.clone(),
                all_limits_passed: false,
                single_limit_reports: Vec::new(),
            };

            let graph_id = self
                .data
                .graph_terms
                .iter()
                .enumerate()
                .find(|(_, term)| term.graph.name == graph_name)
                .ok_or_else(|| eyre!("Graph name '{}' not found in integrand", graph_name))?
                .0;

            for limit in limits {
                let single_limit_report = self.test_single_ir_limit_impl(
                    graph_id,
                    &limit,
                    &mut rng,
                    ir_profile_settings,
                    model,
                )?;

                graph_report.single_limit_reports.push(single_limit_report);
            }

            graph_report.all_limits_passed = graph_report
                .single_limit_reports
                .iter()
                .all(|report| report.passed);

            result.results_per_graph.push(graph_report);
        }

        result.all_passed = result
            .results_per_graph
            .iter()
            .all(|graph_report| graph_report.all_limits_passed);

        Ok(result)
    }

    fn enumerate_ir_limits(&self) -> Vec<(String, Vec<IrLimit>)> {
        self.data
            .graph_terms
            .iter()
            .map(|term| {
                let graph_name = term.graph.name.clone();
                let limits = term.enumerate_ir_limits();
                (graph_name, limits)
            })
            .collect()
    }

    fn parse_select_limits_and_graphs(&self, input: &str) -> Result<Vec<(String, Vec<IrLimit>)>> {
        input
            .split(';')
            .map(|graph_info_string| {
                let mut parts = graph_info_string.split(' ');
                let graph_name = parts
                    .next()
                    .ok_or_else(|| eyre!("Expected graph name in select_limits_and_graphs"))?
                    .to_string();

                let is_valid_name = self
                    .data
                    .graph_terms
                    .iter()
                    .any(|term| term.graph.name == graph_name);

                if !is_valid_name {
                    return Err(eyre!(
                        "Graph name '{}' in select_limits_and_graphs does not match any graph in the integrand",
                        graph_name
                    ));
                }

                let ir_limits = parts
                    .map(IrLimit::parse_limit)
                    .collect::<Result<Vec<_>, _>>()?;

                if ir_limits.is_empty() {
                    return Err(eyre!(
                        "No IR limits specified for graph '{}' in select_limits_and_graphs",
                        graph_name
                    ));
                }

                if ir_limits.iter().any(|limit| !limit.is_valid(self.data.graph_terms.iter().find(|term| term.graph.name == graph_name).unwrap().graph.loop_momentum_basis.loop_edges.len()).is_ok()) {
                    return Err(eyre!(
                        "One or more IR limits specified for graph '{}' in select_limits_and_graphs are not valid",
                        graph_name
                    ));
                }

                Ok((graph_name, ir_limits))
            })
            .collect::<Result<Vec<_>, _>>()
    }

    fn test_single_ir_limit_impl(
        &mut self,
        graph_id: usize,
        ir_limit: &IrLimit,
        rng: &mut MonteCarloRng,
        approach_settings: &IRProfileSetting,
        model: &Model,
    ) -> Result<SingleLimitReport> {
        let edges_in_limit = ir_limit.get_all_edges()?;
        // find lmb
        let lmb = self.data.graph_terms[graph_id].lmb_with_loop_edges(edges_in_limit.as_slice())?;

        // let model_parameter_cache = model.generate_values::<f128>();

        let momenta = ir_limit.get_momenta(rng, &self.settings, approach_settings);
        let non_limit_loops = lmb
            .loop_edges
            .iter_enumerated()
            .filter_map(|(loop_id, edge_id)| {
                if !edges_in_limit.contains(edge_id) {
                    Some(loop_id)
                } else {
                    None
                }
            })
            .collect_vec();

        let non_limit_momenta = non_limit_loops
            .iter()
            .map(|loop_id| (*loop_id, sample_random_unit_vector(rng)))
            .collect_vec();
        let externals = self.data.graph_terms[graph_id]
            .graph
            .get_external_signature();

        let dependent_momenta_constructor = DependentMomentaConstructor::Amplitude(&externals);

        let loop_number = lmb.loop_edges.len();

        let mut limit_data = LimitData { data: Vec::new() };

        for (loop_mom_id, lambda_point) in momenta.into_iter().enumerate() {
            let mut loop_moms: LoopMomenta<F<_>> = (0..loop_number)
                .map(|_| ThreeMomentum::new(F::from_f64(0.0), F::from_f64(0.0), F::from_f64(0.0)))
                .collect();

            for (loop_id, momentum) in non_limit_momenta.iter() {
                loop_moms[*loop_id] = *momentum;
            }

            for tagged_momenta in &lambda_point.momenta {
                let edge_id = tagged_momenta.tag;
                let loop_id = lmb
                    .loop_edges
                    .iter()
                    .position(|loop_edge| loop_edge == &edge_id)
                    .unwrap_or_else(|| {
                        unreachable!("corrupted lmb and ir limit: {}", ir_limit);
                    });

                loop_moms[LoopIndex(loop_id)] = tagged_momenta.momentum;
            }

            let sample_in_cmb = MomentumSample::new(
                loop_moms,
                loop_mom_id,
                &self.settings.kinematics.externals,
                0,
                F::from_f64(1.0),
                dependent_momenta_constructor,
                None,
            )?;

            let sample = sample_in_cmb.lmb_transform(
                &lmb,
                &self.data.graph_terms[graph_id].graph.loop_momentum_basis,
            );

            let sample_flattened = sample
                .loop_moms()
                .iter()
                .flat_map(|mom| [mom.px, mom.py, mom.pz])
                .collect_vec();

            let symbolica_sample = Sample::Discrete(
                F(1.0),
                graph_id,
                Some(Box::new(Sample::Continuous(F(1.0), sample_flattened))),
            );

            limit_data.data.push(LambdaPointEval {
                lambda_point,
                value: self
                    .evaluate_sample(
                        &symbolica_sample,
                        model,
                        sample.one(),
                        0,
                        false,
                        Complex::new_re(F::from_f64(100.0 * self.settings.kinematics.e_cm)),
                    )
                    .unwrap(),
            });
        }

        let slope = limit_data.extract_power()?;

        let num_soft = ir_limit.num_soft();
        let scaling = slope.exponent + ((num_soft * 3) as f64);
        let passed = scaling > 0.0;

        let result = SingleLimitReport {
            limit_name: format!("{}", ir_limit),
            num_soft: ir_limit.num_soft(),
            scaling,
            passed,
            power_law_fit: slope,
        };

        Ok(result)
    }
}

impl CrossSectionIntegrand {
    pub fn ir_profile_completion_entries(&self) -> Vec<(String, Vec<String>)> {
        self.enumerate_ir_limits()
            .into_iter()
            .map(|(graph_name, limits)| {
                (
                    graph_name,
                    limits.into_iter().map(|limit| limit.to_string()).collect(),
                )
            })
            .collect()
    }

    pub fn test_ir(
        &mut self,
        ir_profile_settings: &IRProfileSetting,
        model: &Model,
    ) -> Result<IrLimitTestReport> {
        // override the sampling to be in momentum space
        self.settings.sampling = SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
            sample_orientations: false,
            sampling_type: DiscreteGraphSamplingType::Default(ParameterizationSettings {
                mode: ParameterizationMode::MomentumSpace,
                mapping: ParameterizationMapping::default(),
                b: 10.0,
            }),
        });

        self.warm_up(model)?;

        let mut rng = MonteCarloRng::new(ir_profile_settings.seed, 0);

        let limits_to_check =
            if let Some(select_limits_and_graphs) = &ir_profile_settings.select_limits_and_graphs {
                self.parse_select_limits_and_graphs(select_limits_and_graphs)?
            } else {
                self.enumerate_ir_limits()
            };

        let mut result = IrLimitTestReport {
            all_passed: false,
            results_per_graph: Vec::new(),
        };

        for (graph_name, limits) in limits_to_check {
            let mut graph_report = GraphIRLimitReport {
                graph_name: graph_name.clone(),
                all_limits_passed: false,
                single_limit_reports: Vec::new(),
            };

            let graph_id = self
                .data
                .graph_terms
                .iter()
                .enumerate()
                .find(|(_, term)| term.graph.name == graph_name)
                .ok_or_else(|| eyre!("Graph name '{}' not found in integrand", graph_name))?
                .0;

            for limit in limits {
                let single_limit_report = self.test_single_ir_limit_impl(
                    graph_id,
                    &limit,
                    &mut rng,
                    ir_profile_settings,
                    model,
                )?;

                graph_report.single_limit_reports.push(single_limit_report);
            }

            graph_report.all_limits_passed = graph_report
                .single_limit_reports
                .iter()
                .all(|report| report.passed);

            result.results_per_graph.push(graph_report);
        }

        result.all_passed = result
            .results_per_graph
            .iter()
            .all(|graph_report| graph_report.all_limits_passed);

        Ok(result)
    }

    fn enumerate_ir_limits(&self) -> Vec<(String, Vec<IrLimit>)> {
        self.data
            .graph_terms
            .iter()
            .map(|term| {
                let graph_name = term.graph.name.clone();
                let limits = term.enumerate_ir_limits();
                (graph_name, limits)
            })
            .collect()
    }

    fn parse_select_limits_and_graphs(&self, input: &str) -> Result<Vec<(String, Vec<IrLimit>)>> {
        input
            .split(';')
            .map(|graph_info_string| {
                let mut parts = graph_info_string.split(' ');
                let graph_name = parts
                    .next()
                    .ok_or_else(|| eyre!("Expected graph name in select_limits_and_graphs"))?
                    .to_string();

                let is_valid_name = self
                    .data
                    .graph_terms
                    .iter()
                    .any(|term| term.graph.name == graph_name);

                if !is_valid_name {
                    return Err(eyre!(
                        "Graph name '{}' in select_limits_and_graphs does not match any graph in the integrand",
                        graph_name
                    ));
                }

                let ir_limits = parts
                    .map(IrLimit::parse_limit)
                    .collect::<Result<Vec<_>, _>>()?;

                if ir_limits.is_empty() {
                    return Err(eyre!(
                        "No IR limits specified for graph '{}' in select_limits_and_graphs",
                        graph_name
                    ));
                }

                if ir_limits.iter().any(|limit| !limit.is_valid(self.data.graph_terms.iter().find(|term| term.graph.name == graph_name).unwrap().graph.loop_momentum_basis.loop_edges.len()).is_ok()) {
                    return Err(eyre!(
                        "One or more IR limits specified for graph '{}' in select_limits_and_graphs are not valid",
                        graph_name
                    ));
                }

                Ok((graph_name, ir_limits))
            })
            .collect::<Result<Vec<_>, _>>()
    }

    fn test_single_ir_limit_impl(
        &mut self,
        graph_id: usize,
        ir_limit: &IrLimit,
        rng: &mut MonteCarloRng,
        approach_settings: &IRProfileSetting,
        model: &Model,
    ) -> Result<SingleLimitReport> {
        let edges_in_limit = ir_limit.get_all_edges()?;

        // find cut that for that as all edges of the limit
        let (cut_id, _esurface) = self.data.graph_terms[graph_id]
            .cut_esurface
            .iter_enumerated()
            .find(|(_cut_id, esurface)| {
                edges_in_limit
                    .iter()
                    .all(|edge| esurface.energies.contains(edge))
            })
            .ok_or(eyre!(
                "could not find cut with all edges of the limit: {}",
                ir_limit
            ))?;

        let cs_cut = &self.data.graph_terms[graph_id].cuts[cut_id];

        let edges_to_flip = cs_cut
            .cut
            .iter_edges(&self.data.graph_terms[graph_id].graph.underlying)
            .map(|(or, _)| or)
            .zip(
                self.data.graph_terms[graph_id]
                    .graph
                    .underlying
                    .iter_edges_of(&cs_cut.cut)
                    .map(|x| x.1),
            )
            .filter_map(|(orientation, edge_id)| {
                if edges_in_limit.contains(&edge_id) && matches!(orientation, Orientation::Reversed)
                {
                    Some(edge_id)
                } else {
                    None
                }
            })
            .collect_vec();

        // find lmb
        let lmb = self.data.graph_terms[graph_id].lmb_with_loop_edges(edges_in_limit.as_slice())?;

        // let model_parameter_cache = model.generate_values::<f128>();

        let momenta = ir_limit.get_momenta(rng, &self.settings, approach_settings);
        let non_limit_loops = lmb
            .loop_edges
            .iter_enumerated()
            .filter_map(|(loop_id, edge_id)| {
                if !edges_in_limit.contains(edge_id) {
                    Some(loop_id)
                } else {
                    None
                }
            })
            .collect_vec();

        let non_limit_momenta = non_limit_loops
            .iter()
            .map(|loop_id| (*loop_id, sample_random_unit_vector(rng)))
            .collect_vec();

        let dependent_momenta_constructor = DependentMomentaConstructor::CrossSection;

        let loop_number = lmb.loop_edges.len();

        let mut limit_data = LimitData { data: Vec::new() };

        for (loop_mom_id, lambda_point) in momenta.into_iter().enumerate() {
            let mut loop_moms: LoopMomenta<F<_>> = (0..loop_number)
                .map(|_| ThreeMomentum::new(F::from_f64(0.0), F::from_f64(0.0), F::from_f64(0.0)))
                .collect();

            for (loop_id, momentum) in non_limit_momenta.iter() {
                loop_moms[*loop_id] = *momentum;
            }

            for tagged_momenta in &lambda_point.momenta {
                let edge_id = tagged_momenta.tag;
                let loop_id = lmb
                    .loop_edges
                    .iter()
                    .position(|loop_edge| loop_edge == &edge_id)
                    .unwrap_or_else(|| {
                        unreachable!("corrupted lmb and ir limit: {}", ir_limit);
                    });

                if edges_to_flip.contains(&edge_id) {
                    loop_moms[LoopIndex(loop_id)] = -tagged_momenta.momentum;
                } else {
                    loop_moms[LoopIndex(loop_id)] = tagged_momenta.momentum
                };
            }

            let sample_in_cmb = MomentumSample::new(
                loop_moms,
                loop_mom_id,
                &self.settings.kinematics.externals,
                0,
                F::from_f64(1.0),
                dependent_momenta_constructor,
                None,
            )?;

            let sample = sample_in_cmb.lmb_transform(
                &lmb,
                &self.data.graph_terms[graph_id].graph.loop_momentum_basis,
            );

            let sample_flattened = sample
                .loop_moms()
                .iter()
                .flat_map(|mom| [mom.px, mom.py, mom.pz])
                .collect_vec();

            let symbolica_sample = Sample::Discrete(
                F(1.0),
                graph_id,
                Some(Box::new(Sample::Continuous(F(1.0), sample_flattened))),
            );

            limit_data.data.push(LambdaPointEval {
                lambda_point,
                value: self
                    .evaluate_sample(
                        &symbolica_sample,
                        model,
                        sample.one(),
                        0,
                        false,
                        Complex::new_re(F::from_f64(100.0 * self.settings.kinematics.e_cm)),
                    )
                    .unwrap(),
            });
        }

        let slope = limit_data.extract_power()?;

        let num_soft = ir_limit.num_soft();
        let scaling = slope.exponent + ((num_soft * 3) as f64);
        let passed = scaling > 0.0;

        let result = SingleLimitReport {
            limit_name: format!("{}", ir_limit),
            num_soft: ir_limit.num_soft(),
            scaling,
            passed,
            power_law_fit: slope,
        };

        Ok(result)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct IrLimit {
    colinear: Vec<Vec<HardOrSoft>>,
    soft: Vec<EdgeIndex>,
}

enum MomentumBuilder<T: FloatLike> {
    Colinear {
        edge_id: EdgeIndex,
        x: F<T>,
        colinear_direction: ThreeMomentum<F<T>>,
        perpendicular_direction: ThreeMomentum<F<T>>,
        is_soft: bool,
    },
    Soft {
        edge_id: EdgeIndex,
        direction: ThreeMomentum<F<T>>,
    },
}

impl Display for IrLimit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for colinear_set in &self.colinear {
            write!(f, "C[")?;

            let mut iter = colinear_set.iter();
            if let Some(last) = iter.next_back() {
                for item in iter {
                    write!(f, "{},", item)?;
                }
                write!(f, "{}]", last)?;
            } else {
                write!(f, "]")?;
            }
        }

        for soft in self.soft.iter() {
            write!(f, "S({})", soft)?;
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum HardOrSoft {
    Hard(EdgeIndex),
    Soft(EdgeIndex),
}

impl HardOrSoft {
    fn index(&self) -> EdgeIndex {
        match self {
            HardOrSoft::Hard(index) => *index,
            HardOrSoft::Soft(index) => *index,
        }
    }
}

impl Display for HardOrSoft {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HardOrSoft::Hard(index) => write!(f, "{}", index),
            HardOrSoft::Soft(index) => write!(f, "S({})", index),
        }
    }
}

impl IrLimit {
    fn canonize(&mut self) {
        for colinear_set in &mut self.colinear.iter_mut() {
            colinear_set.sort();
        }

        self.colinear.sort();
        self.soft.sort();
    }

    fn new_pure_colinear(colinear_edges: Vec<EdgeIndex>) -> Self {
        let colinear = vec![colinear_edges.into_iter().map(HardOrSoft::Hard).collect()];

        let mut result = IrLimit {
            colinear,
            soft: Vec::new(),
        };
        result.canonize();
        result
    }

    fn new_pure_soft(soft_edges: Vec<EdgeIndex>) -> Self {
        let mut result = IrLimit {
            colinear: Vec::new(),
            soft: soft_edges,
        };
        result.canonize();
        result
    }

    fn num_soft(&self) -> usize {
        self.colinear
            .iter()
            .flatten()
            .filter(|edge| matches!(edge, HardOrSoft::Soft(_)))
            .count()
            + self.soft.len()
    }

    fn check_min_colinear_size(&self) -> bool {
        self.colinear
            .iter()
            .all(|colinear_set| colinear_set.len() >= 2)
    }

    fn is_valid(&self, loop_number: usize) -> Result<()> {
        if !self.check_min_colinear_size() {
            return Err(eyre!("colinear sets must have at least two edges"));
        }

        let all_edges = self.get_all_edges()?;

        if all_edges.len() > loop_number {
            return Err(eyre!("not enough degrees of freedom to setup IR limit"));
        }

        Ok(())
    }

    fn get_all_edges(&self) -> Result<Vec<EdgeIndex>> {
        let colinear_edges = self
            .colinear
            .iter()
            .flatten()
            .map(HardOrSoft::index)
            .collect_vec();
        let soft_edges = self.soft.iter().copied().collect_vec();

        let all_edges: Vec<EdgeIndex> = colinear_edges
            .into_iter()
            .chain(soft_edges)
            .sorted()
            .collect();

        // check for duplicates
        let mut unique_edges = all_edges.clone();
        unique_edges.dedup();

        if unique_edges.len() != all_edges.len() {
            return Err(eyre!("Edges specified in ir limit must be unique")); // duplicates found
        }

        Ok(all_edges)
    }

    fn parse_limit(limit: &str) -> Result<IrLimit> {
        let mut colinear_sets = Vec::new();
        let mut soft_edges = Vec::new();

        let mut char_iter = limit.chars().enumerate();

        while let Some((char_position, char)) = char_iter.next() {
            match char {
                'C' => {
                    let (_opening_bracket_position, opening_bracket) =
                        char_iter.next().ok_or_else(|| {
                            eyre!(
                                "Expected opening bracket after 'C' at position {}",
                                char_position
                            )
                        })?;

                    if opening_bracket != '[' {
                        return Err(eyre!(
                            "Expected '[' after 'C' at position {}, found '{}'",
                            char_position,
                            opening_bracket
                        ));
                    }

                    let mut colinear_set_str = String::new();

                    let mut closing_bracket_found = false;
                    for (_next_char_position, next_char) in char_iter.by_ref() {
                        if next_char == ']' {
                            closing_bracket_found = true;
                            break;
                        }
                        colinear_set_str.push(next_char);
                    }

                    if !closing_bracket_found {
                        return Err(eyre!(
                            "Expected closing bracket ']' for colinear set at position {}",
                            char_position
                        ));
                    }

                    let edges = colinear_set_str.trim().split(',');

                    let mut colinear_set = Vec::new();

                    for edge in edges {
                        let trimmed_edge = edge.trim();
                        if trimmed_edge.is_empty() {
                            return Err(eyre!(
                                "Empty edge found in colinear set at position {}",
                                char_position
                            ));
                        }

                        if trimmed_edge.starts_with('S') {
                            let mut trimmed_edge_iter = trimmed_edge.chars().skip(1);
                            let opening_bracket = trimmed_edge_iter
                                .next()
                                .ok_or(eyre!("Expected '(' after 'S' in soft edge at position"))?;

                            if opening_bracket != '(' {
                                return Err(eyre!(
                                    "Expected '(' after 'S' in soft edge at position , found ''",
                                ));
                            }

                            let mut edge_str = String::new();
                            let mut closing_bracket_found = false;

                            for next_char in trimmed_edge_iter {
                                if next_char == ')' {
                                    closing_bracket_found = true;
                                    break;
                                }
                                edge_str.push(next_char);
                            }

                            if !closing_bracket_found {
                                return Err(eyre!(
                                    "Expected closing bracket ')' for soft edge at position {}",
                                    char_position
                                ));
                            }

                            let edge_index = Self::parse_edge(&edge_str)?;
                            colinear_set.push(HardOrSoft::Soft(edge_index));
                        } else {
                            let edge_index = Self::parse_edge(trimmed_edge)?;
                            colinear_set.push(HardOrSoft::Hard(edge_index));
                        }
                    }
                    colinear_sets.push(colinear_set);
                }
                'S' => {
                    let (_opening_bracket_position, opening_bracket) =
                        char_iter.next().ok_or_else(|| {
                            eyre!(
                                "Expected opening bracket '(' after 'S' at position {}",
                                char_position
                            )
                        })?;

                    if opening_bracket != '(' {
                        return Err(eyre!(
                            "Expected '(' after 'S' at position {}, found '{}'",
                            char_position,
                            opening_bracket
                        ));
                    }

                    let mut edge_str = String::new();
                    let mut closing_bracket_found = false;

                    for (_next_char_position, next_char) in char_iter.by_ref() {
                        if next_char == ')' {
                            closing_bracket_found = true;
                            break;
                        }
                        edge_str.push(next_char);
                    }

                    if !closing_bracket_found {
                        return Err(eyre!(
                            "Expected closing bracket ')' for soft edge at position {}",
                            char_position
                        ));
                    }

                    let edge_index = Self::parse_edge(&edge_str)?;
                    soft_edges.push(edge_index);
                }
                _ => {
                    return Err(eyre!(
                        "Unexpected character '{}' at position {}",
                        char,
                        char_position
                    ));
                }
            }
        }

        let mut ir_limit = IrLimit {
            colinear: colinear_sets,
            soft: soft_edges,
        };

        ir_limit.canonize();

        Ok(ir_limit)
    }

    fn parse_edge(edge: &str) -> Result<EdgeIndex> {
        let mut edge = String::from(edge);
        edge = String::from(edge.trim());

        if edge.len() < 2 {
            return Err(eyre!("Edge must be at least two characters long"));
        }

        if edge.remove(0) != 'e' {
            return Err(eyre!("Edge must start with 'e'"));
        }

        let edge_id: usize = edge
            .parse()
            .map_err(|_| eyre!("Edge must be a valid integer, got: {}", edge))?;

        Ok(EdgeIndex::from(edge_id))
    }

    fn get_momentum_builders(&self, rng: &mut MonteCarloRng) -> Vec<MomentumBuilder<f64>> {
        let mut momentum_builder = Vec::new();

        for colinear_set in &self.colinear {
            let direction_for_set: ThreeMomentum<F<f64>> = sample_random_unit_vector(rng);

            let x_variables: Vec<F<f64>> = (0..colinear_set.len())
                .map(|_| F::from_f64(rng.random::<f64>() * 0.8 + 0.1))
                .sorted_by(|a, b| a.partial_cmp(b).unwrap())
                .collect_vec();

            for (edge, x) in colinear_set.iter().zip(x_variables) {
                let edge_id = edge.index();
                let direction: ThreeMomentum<F<f64>> = sample_random_unit_vector(rng);

                let perpendicular = direction - direction * (direction * direction_for_set);

                let perpendicular_norm = perpendicular.norm();
                let perpendicular = perpendicular * perpendicular_norm.inv();

                let is_soft = matches!(edge, HardOrSoft::Soft(_));

                momentum_builder.push(MomentumBuilder::Colinear {
                    edge_id,
                    x,
                    colinear_direction: direction_for_set,
                    perpendicular_direction: perpendicular,
                    is_soft,
                });
            }
        }

        for soft_edge in &self.soft {
            let direction = sample_random_unit_vector(rng);
            momentum_builder.push(MomentumBuilder::Soft {
                edge_id: *soft_edge,
                direction,
            });
        }

        momentum_builder
    }

    fn get_momenta(
        &self,
        rng: &mut MonteCarloRng,
        settings: &RuntimeSettings,
        approach_settings: &IRProfileSetting,
    ) -> Vec<LambdaPoint<f64>> {
        let momentum_builders = self.get_momentum_builders(rng);

        let lambda_values: Vec<F<f64>> = logspace(
            approach_settings.lambda_exp_start,
            approach_settings.lambda_exp_end,
            approach_settings.steps,
            10.0,
        )
        .into_iter()
        .map(F::from_f64)
        .collect();

        lambda_values
            .into_iter()
            .map(|lambda| LambdaPoint {
                momenta: momentum_builders
                    .iter()
                    .map(|builder| match builder {
                        MomentumBuilder::Colinear {
                            edge_id,
                            x,
                            colinear_direction,
                            perpendicular_direction,
                            is_soft,
                        } => {
                            let momentum = if *is_soft {
                                (colinear_direction * x + perpendicular_direction * lambda)
                                    * F::from_f64(settings.kinematics.e_cm)
                                    * lambda
                            } else {
                                (colinear_direction * x + perpendicular_direction * lambda)
                                    * F::from_f64(settings.kinematics.e_cm)
                            };
                            TaggedMomenta {
                                momentum,
                                tag: *edge_id,
                            }
                        }
                        MomentumBuilder::Soft { edge_id, direction } => {
                            let momentum =
                                direction * lambda * F::from_f64(settings.kinematics.e_cm);
                            TaggedMomenta {
                                momentum,
                                tag: *edge_id,
                            }
                        }
                    })
                    .collect(),
                lambda,
            })
            .collect()
    }
}

fn sample_random_unit_vector<T: FloatLike>(rng: &mut MonteCarloRng) -> ThreeMomentum<F<T>> {
    let x_1 = F::<T>::from_f64(rng.random::<f64>());
    let x_2 = F::<T>::from_f64(rng.random::<f64>());
    let x_3 = F::<T>::from_f64(rng.random::<f64>());
    let x_4 = F::<T>::from_f64(rng.random::<f64>());

    let (k_x, k_y) = box_muller(x_1, x_2);
    let k_z = box_muller(x_3, x_4).0;

    let unnormalized_momentum = ThreeMomentum::new(k_x, k_y, k_z);
    let norm = unnormalized_momentum.norm();
    unnormalized_momentum * norm.inv()
}

struct TaggedMomenta<T> {
    momentum: ThreeMomentum<T>,
    tag: EdgeIndex,
}

struct LambdaPoint<T: FloatLike> {
    lambda: F<T>,
    momenta: Vec<TaggedMomenta<F<T>>>,
}

struct LambdaPointEval<T: FloatLike> {
    lambda_point: LambdaPoint<T>,
    value: EvaluationResult,
}

struct LimitData<T: FloatLike> {
    data: Vec<LambdaPointEval<T>>,
}

impl<T: FloatLike> LimitData<T> {
    fn extract_power(&self) -> Result<PowerLawFit> {
        let x = self
            .data
            .iter()
            .map(|point_eval| point_eval.lambda_point.lambda.to_f64())
            .collect_vec();

        let y = self
            .data
            .iter()
            .map(|point_eval| point_eval.value.integrand_result.re.to_f64())
            .collect_vec();

        let result = fit_power_law(x.clone(), y.clone())?;

        if result.r_squared < 0.9 {
            warn!("low r^2 value found for input data");
            warn!("x: {:?}", x);
            warn!("y: {:?}", y);
        }

        Ok(result)
    }
}

#[derive(Debug, Clone)]
pub struct PowerLawFit {
    exponent: f64,
    coefficient: f64,
    constant_offset: f64,
    r_squared: f64,
}

// using varpro fit  y =  a x^p + c
fn fit_power_law(x: Vec<f64>, y: Vec<f64>) -> Result<PowerLawFit> {
    if x.len() != y.len() {
        return Err(eyre!(
            "fit_power_law requires x and y to have the same length"
        ));
    }
    if x.len() < 2 {
        return Err(eyre!("fit_power_law requires at least two observations"));
    }
    if x.iter().any(|value| !value.is_finite() || *value <= 0.0) {
        return Err(eyre!(
            "fit_power_law requires strictly positive, finite x values"
        ));
    }
    if y.iter().any(|value| !value.is_finite()) {
        return Err(eyre!("fit_power_law requires finite y values"));
    }

    let x_scale = (x.iter().map(|value| value.ln()).sum::<f64>() / x.len() as f64).exp();
    if !x_scale.is_finite() || x_scale <= 0.0 {
        return Err(eyre!(
            "fit_power_law could not determine a valid x rescaling factor"
        ));
    }

    let x_normalized = x.iter().map(|value| value / x_scale).collect::<Vec<_>>();

    let x_data = DVector::from_vec(x_normalized.clone());
    let y_data = DVector::from_vec(y.clone());

    let initial_exponent = {
        let mut log_x = Vec::new();
        let mut log_y_shifted = Vec::new();
        let min_y = y.iter().copied().fold(f64::INFINITY, f64::min);
        let y_shift = if min_y <= 0.0 { 1.0 - min_y } else { 0.0 };
        for (xv, yv) in x.iter().zip(y.iter()) {
            let shifted = yv + y_shift;
            if shifted > 0.0 {
                log_x.push(xv.ln());
                log_y_shifted.push(shifted.ln());
            }
        }

        if log_x.len() >= 2 {
            let mean_x = log_x.iter().sum::<f64>() / log_x.len() as f64;
            let mean_y = log_y_shifted.iter().sum::<f64>() / log_y_shifted.len() as f64;
            let covariance = log_x
                .iter()
                .zip(log_y_shifted.iter())
                .map(|(xv, yv)| (xv - mean_x) * (yv - mean_y))
                .sum::<f64>();
            let variance = log_x.iter().map(|xv| (xv - mean_x).powi(2)).sum::<f64>();
            if variance > 0.0 {
                covariance / variance
            } else {
                -1.0
            }
        } else {
            -1.0
        }
    };

    const LARGE_FINITE_GUARD: f64 = 1.0e200;

    let model = SeparableModelBuilder::new(["p"])
        .initial_parameters(vec![initial_exponent])
        .independent_variable(x_data)
        .function(["p"], |x: &DVector<f64>, p: f64| {
            x.map(|xv| {
                let value = xv.powf(p);
                if value.is_finite() {
                    value
                } else if value.is_nan() {
                    LARGE_FINITE_GUARD
                } else if value.is_sign_negative() {
                    -LARGE_FINITE_GUARD
                } else {
                    LARGE_FINITE_GUARD
                }
            })
        })
        .partial_deriv("p", |x: &DVector<f64>, p: f64| {
            x.map(|xv| {
                let deriv = xv.powf(p) * xv.ln();
                if deriv.is_finite() {
                    deriv
                } else if deriv.is_nan() {
                    LARGE_FINITE_GUARD
                } else if deriv.is_sign_negative() {
                    -LARGE_FINITE_GUARD
                } else {
                    LARGE_FINITE_GUARD
                }
            })
        })
        .invariant_function(|x: &DVector<f64>| DVector::from_element(x.len(), 1.0))
        .build()
        .map_err(|error| eyre!("could not build varpro model for power-law fit: {error}"))?;

    let problem = SeparableProblemBuilder::new(model)
        .observations(y_data)
        .build()
        .map_err(|error| eyre!("could not build varpro problem for power-law fit: {error}"))?;

    let fit_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        LevMarSolver::default().solve(problem)
    }))
    .map_err(|_| eyre!("varpro/nalgebra panicked while solving power-law fit"))?
    .map_err(|error| eyre!("varpro solve failed for power-law fit: {:?}", error))?;

    if !fit_result.minimization_report.termination.was_successful() {
        return Err(eyre!(
            "varpro did not terminate successfully in power-law fit"
        ));
    }

    let exponent = fit_result.nonlinear_parameters()[0];
    let coefficients = fit_result
        .linear_coefficients()
        .ok_or_else(|| eyre!("varpro did not return linear coefficients"))?;
    let coefficient_normalized = coefficients[0];
    let constant_offset = coefficients[1];

    let exponent = if exponent.is_finite() {
        exponent
    } else {
        return Err(eyre!("power-law fit returned non-finite exponent"));
    };

    let coefficient_scale = x_scale.powf(-exponent);
    if !coefficient_scale.is_finite() {
        return Err(eyre!(
            "power-law fit produced non-finite coefficient rescaling"
        ));
    }

    let coefficient = coefficient_normalized * coefficient_scale;
    if !coefficient.is_finite() || !constant_offset.is_finite() {
        return Err(eyre!(
            "power-law fit returned non-finite linear coefficients"
        ));
    }

    let mean_y = y.iter().sum::<f64>() / y.len() as f64;
    let ss_tot = y.iter().map(|yv| (yv - mean_y).powi(2)).sum::<f64>();
    let ss_res = x
        .iter()
        .zip(y.iter())
        .map(|(xv, yv)| {
            let prediction = coefficient * xv.powf(&exponent) + constant_offset;
            (yv - prediction).powi(2)
        })
        .sum::<f64>();
    let r_squared = if ss_tot > 0.0 {
        1.0 - ss_res / ss_tot
    } else {
        1.0
    };

    Ok(PowerLawFit {
        exponent,
        coefficient,
        constant_offset,
        r_squared,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_display() {
        let ir_limit = IrLimit {
            colinear: vec![
                vec![
                    HardOrSoft::Hard(EdgeIndex::from(1)),
                    HardOrSoft::Hard(EdgeIndex::from(2)),
                    HardOrSoft::Hard(EdgeIndex::from(3)),
                ],
                vec![
                    HardOrSoft::Hard(EdgeIndex::from(4)),
                    HardOrSoft::Soft(EdgeIndex::from(5)),
                ],
            ],
            soft: vec![EdgeIndex::from(6), EdgeIndex::from(7)],
        };

        let display = ir_limit.to_string();
        let expected = "C[e1,e2,e3]C[e4,S(e5)]S(e6)S(e7)";

        assert_eq!(display, expected);
    }

    #[test]
    fn parse_edge() {
        let edge_str = "e5";
        let edge_index = IrLimit::parse_edge(edge_str).unwrap();
        assert_eq!(edge_index, EdgeIndex::from(5));

        let invalid_edge_str = "5"; // missing 'e'
        assert!(IrLimit::parse_edge(invalid_edge_str).is_err());

        let invalid_edge_str2 = "e"; // too short
        assert!(IrLimit::parse_edge(invalid_edge_str2).is_err());

        let invalid_edge_str3 = "e5a"; // not a valid integer
        assert!(IrLimit::parse_edge(invalid_edge_str3).is_err());
    }

    #[test]
    fn parse_limit() {
        let limit_str = "C[e1,e2,e3]C[e4,S(e5)]S(e6)S(e7)";
        let ir_limit = IrLimit::parse_limit(limit_str).unwrap();

        assert_eq!(ir_limit.colinear.len(), 2, "Expected two colinear sets");
        assert_eq!(ir_limit.soft.len(), 2, "Expected two soft edges");

        assert_eq!(
            ir_limit.colinear[0],
            vec![
                HardOrSoft::Hard(EdgeIndex::from(1)),
                HardOrSoft::Hard(EdgeIndex::from(2)),
                HardOrSoft::Hard(EdgeIndex::from(3))
            ],
            "First colinear set does not match"
        );
        assert_eq!(
            ir_limit.colinear[1],
            vec![
                HardOrSoft::Hard(EdgeIndex::from(4)),
                HardOrSoft::Soft(EdgeIndex::from(5))
            ],
            "Second colinear set does not match"
        );

        assert_eq!(
            ir_limit.soft,
            vec![EdgeIndex::from(6), EdgeIndex::from(7)],
            "Soft edges do not match"
        );

        let invalid_limit_str = "C[e1,e2,e3]C[e4,S(e5)]S(e6)S(e7, e8)";
        assert!(
            IrLimit::parse_limit(invalid_limit_str).is_err(),
            "Expected error"
        );

        let invalid_limit_str2 = "C[e1,e2,e3C[e4,S(e5)]S(e6)";
        assert!(
            IrLimit::parse_limit(invalid_limit_str2).is_err(),
            "Expected error for unmatched brackets"
        );
    }

    #[test]
    fn fit_power_law_recovers_known_parameters() {
        let exponent = -1.75;
        let coefficient = 3.2;
        let offset = 0.6;

        let x = vec![0.2, 0.35, 0.5, 0.8, 1.1, 1.7, 2.4, 3.3];
        let y = x
            .iter()
            .map(|xv| coefficient * xv.powf(&exponent) + offset)
            .collect::<Vec<_>>();

        let fit = fit_power_law(x, y).expect("power-law fit should succeed");

        assert!((fit.exponent - exponent).abs() < 1e-3);
        assert!((fit.coefficient - coefficient).abs() < 1e-3);
        assert!(fit.r_squared > 0.999_999);
    }

    #[test]
    fn fit_power_law_on_panicking_data() {
        let x = vec![
            0.001,
            0.0007847599703514615,
            0.0006158482110660267,
            0.0004832930238571752,
            0.000379269019073225,
            0.00029763514416313193,
            0.00023357214690901214,
            0.00018329807108324357,
            0.0001438449888287663,
            0.00011288378916846895,
            8.858667904100833e-5,
            6.951927961775606e-5,
            5.4555947811685143e-5,
            4.281332398719396e-5,
            3.359818286283781e-5,
            2.6366508987303556e-5,
            2.06913808111479e-5,
            1.6237767391887242e-5,
            1.2742749857031348e-5,
            1e-5,
        ];

        let y = vec![
            0.002950535397303611,
            0.004802153664059006,
            0.007811787818354787,
            0.012702615924354177,
            0.020646917328122072,
            0.033559978182893246,
            0.05451887019444257,
            0.08857211912982166,
            0.14386785496026278,
            0.2341369427740574,
            0.378563791513443,
            0.6228243708610535,
            1.0381001830101013,
            1.609904408454895,
            2.4859671592712402,
            4.144830703735352,
            7.431371688842773,
            8.509048461914063,
            -28.14263916015625,
            90.61788940429688,
        ];

        let _ = fit_power_law(x, y).is_ok();
    }
}
