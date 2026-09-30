use std::time::Duration;

use serde::{Deserialize, Serialize};

use crate::cff::CffEnergyDegreeBoundReport;

#[derive(Debug, Clone, Copy, Default)]
pub struct EvaluatorBuildTimings {
    pub spenso_time: Duration,
    pub symbolica_time: Duration,
}

impl std::ops::AddAssign for EvaluatorBuildTimings {
    fn add_assign(&mut self, other: Self) {
        self.spenso_time += other.spenso_time;
        self.symbolica_time += other.symbolica_time;
    }
}

/// Durations sum work within graph jobs; they are not process wall-clock time.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenerationTimings {
    #[serde(default)]
    pub evaluator_count: usize,
    #[serde(default)]
    pub total_time: Duration,
    #[serde(default)]
    pub evaluator_spenso_time: Duration,
    #[serde(default, alias = "evaluator_build_time")]
    pub evaluator_symbolica_time: Duration,
    #[serde(default)]
    pub evaluator_compile_time: Duration,
}

impl GenerationTimings {
    pub fn evaluator_build_time(&self) -> Duration {
        self.evaluator_spenso_time + self.evaluator_symbolica_time
    }

    pub fn expression_build_time(&self) -> Duration {
        self.total_time
            .saturating_sub(self.evaluator_build_time())
            .saturating_sub(self.evaluator_compile_time)
    }

    pub fn add_evaluator_build_timings(&mut self, timings: EvaluatorBuildTimings) {
        self.evaluator_spenso_time += timings.spenso_time;
        self.evaluator_symbolica_time += timings.symbolica_time;
    }

    pub fn merge_in_place(&mut self, other: &Self) {
        self.evaluator_count += other.evaluator_count;
        self.total_time += other.total_time;
        self.evaluator_spenso_time += other.evaluator_spenso_time;
        self.evaluator_symbolica_time += other.evaluator_symbolica_time;
        self.evaluator_compile_time += other.evaluator_compile_time;
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RepresentationGenerationStats {
    pub representation: three_dimensional_reps::RepresentationMode,
    #[serde(flatten)]
    pub timings: GenerationTimings,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GraphGenerationStats {
    #[serde(flatten)]
    pub timings: GenerationTimings,
    /// Measured disjoint work, in the requested generation order. Shared graph,
    /// Taylor, integration and geometry work belongs only to the aggregate.
    #[serde(default)]
    pub representations: Vec<RepresentationGenerationStats>,
    /// Generation-time CFF diagnostics intentionally omitted from persisted
    /// generation summaries.
    #[serde(skip)]
    pub cff_energy_degree_bound_reports: Vec<CffEnergyDegreeBoundReport>,
}

impl GraphGenerationStats {
    /// Attribute one completed evaluator build; the enclosing graph timer owns
    /// aggregate wall duration, so only the representation's total is added here.
    pub(crate) fn record_evaluator_build(
        &mut self,
        representation: three_dimensional_reps::RepresentationMode,
        timings: EvaluatorBuildTimings,
        evaluator_count: usize,
        elapsed: Duration,
    ) {
        self.timings.add_evaluator_build_timings(timings);
        self.timings.evaluator_count += evaluator_count;
        let entry = self.representation_mut(representation);
        entry.add_evaluator_build_timings(timings);
        entry.evaluator_count += evaluator_count;
        entry.total_time += elapsed;
    }

    pub fn representation_mut(
        &mut self,
        representation: three_dimensional_reps::RepresentationMode,
    ) -> &mut GenerationTimings {
        let index = self
            .representations
            .iter()
            .position(|stats| stats.representation == representation)
            .unwrap_or_else(|| {
                self.representations.push(RepresentationGenerationStats {
                    representation,
                    timings: GenerationTimings::default(),
                });
                self.representations.len() - 1
            });
        &mut self.representations[index].timings
    }

    /// The aggregate includes all measured work. Representation rows are
    /// subsets, so their complement includes common preparation and overhead.
    pub fn shared_stats(&self) -> GenerationTimings {
        let mut shared = self.timings.clone();
        for representation in &self.representations {
            let timing = &representation.timings;
            shared.evaluator_count = shared
                .evaluator_count
                .saturating_sub(timing.evaluator_count);
            shared.total_time = shared.total_time.saturating_sub(timing.total_time);
            shared.evaluator_spenso_time = shared
                .evaluator_spenso_time
                .saturating_sub(timing.evaluator_spenso_time);
            shared.evaluator_symbolica_time = shared
                .evaluator_symbolica_time
                .saturating_sub(timing.evaluator_symbolica_time);
            shared.evaluator_compile_time = shared
                .evaluator_compile_time
                .saturating_sub(timing.evaluator_compile_time);
        }
        shared
    }

    pub fn evaluator_build_time(&self) -> Duration {
        self.timings.evaluator_build_time()
    }

    pub fn expression_build_time(&self) -> Duration {
        self.timings.expression_build_time()
    }

    pub fn merge_in_place(&mut self, other: &Self) {
        self.timings.merge_in_place(&other.timings);
        for representation in &other.representations {
            self.representation_mut(representation.representation)
                .merge_in_place(&representation.timings);
        }
        for report in &other.cff_energy_degree_bound_reports {
            if !self.cff_energy_degree_bound_reports.contains(report) {
                self.cff_energy_degree_bound_reports.push(report.clone());
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NamedGraphGenerationReport {
    pub integrand_name: String,
    pub graph_name: String,
    pub stats: GraphGenerationStats,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct GeneratedGraphKey {
    pub process_id: usize,
    pub integrand_name: String,
    pub graph_name: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneratedGraphReport {
    pub process_id: usize,
    pub integrand_name: String,
    pub graph_name: String,
    pub stats: GraphGenerationStats,
}

impl GeneratedGraphReport {
    pub fn key(&self) -> GeneratedGraphKey {
        GeneratedGraphKey {
            process_id: self.process_id,
            integrand_name: self.integrand_name.clone(),
            graph_name: self.graph_name.clone(),
        }
    }

    pub fn merge_in_place(&mut self, other: &Self) {
        debug_assert_eq!(self.key(), other.key());
        self.stats.merge_in_place(&other.stats);
    }
}

pub fn merge_generated_graph_reports(
    reports: &mut Vec<GeneratedGraphReport>,
    updates: Vec<GeneratedGraphReport>,
) {
    for update in updates {
        if let Some(existing) = reports
            .iter_mut()
            .find(|report| report.key() == update.key())
        {
            existing.merge_in_place(&update);
        } else {
            reports.push(update);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cff::CffEnergyBoundSourceKind;
    use three_dimensional_reps::RepresentationMode;

    #[test]
    fn representation_timings_preserve_order_shared_costs_and_roundtrip() {
        let mut stats = GraphGenerationStats {
            timings: GenerationTimings {
                total_time: Duration::from_secs(100),
                evaluator_count: 1,
                evaluator_spenso_time: Duration::from_secs(3),
                evaluator_symbolica_time: Duration::from_secs(4),
                ..Default::default()
            },
            ..Default::default()
        };
        for (representation, elapsed) in
            [(RepresentationMode::Ltd, 20), (RepresentationMode::Cff, 30)]
        {
            stats.record_evaluator_build(
                representation,
                EvaluatorBuildTimings {
                    spenso_time: Duration::from_secs(5),
                    symbolica_time: Duration::from_secs(7),
                },
                2,
                Duration::from_secs(elapsed),
            );
        }
        // A later compiler reports modes in another traversal order. Merging
        // must preserve generation order and must not duplicate shared work.
        let mut compilation = GraphGenerationStats::default();
        compilation.timings.total_time = Duration::from_secs(13);
        compilation.timings.evaluator_compile_time = Duration::from_secs(13);
        for (representation, elapsed) in
            [(RepresentationMode::Cff, 5), (RepresentationMode::Ltd, 7)]
        {
            let timing = compilation.representation_mut(representation);
            timing.total_time = Duration::from_secs(elapsed);
            timing.evaluator_compile_time = Duration::from_secs(elapsed);
        }
        stats.merge_in_place(&compilation);
        assert_eq!(
            stats
                .representations
                .iter()
                .map(|row| row.representation)
                .collect::<Vec<_>>(),
            [RepresentationMode::Ltd, RepresentationMode::Cff]
        );
        let shared = stats.shared_stats();
        assert_eq!(shared.total_time, Duration::from_secs(51));
        assert_eq!(shared.evaluator_count, 1);
        assert_eq!(shared.evaluator_spenso_time, Duration::from_secs(3));
        assert_eq!(shared.evaluator_symbolica_time, Duration::from_secs(4));
        assert_eq!(shared.evaluator_compile_time, Duration::from_secs(1));
        let mut reconstructed = shared;
        for row in &stats.representations {
            reconstructed.merge_in_place(&row.timings);
        }
        assert_eq!(reconstructed, stats.timings);
        assert_eq!(
            stats.representations[0].timings.expression_build_time(),
            Duration::from_secs(8)
        );
        let encoded = serde_json::to_value(&stats).unwrap();
        assert!(encoded.get("total_time").is_some());
        assert!(encoded.get("timings").is_none());
        assert_eq!(encoded["representations"][0]["representation"], "ltd");
        let restored: GraphGenerationStats = serde_json::from_value(encoded).unwrap();
        assert_eq!(restored.timings, stats.timings);
        assert_eq!(restored.shared_stats(), stats.shared_stats());
        assert_eq!(
            restored.representations[0].timings,
            stats.representations[0].timings
        );
        assert_eq!(
            restored.representations[1].timings,
            stats.representations[1].timings
        );
    }

    #[test]
    fn cff_energy_bound_reports_are_merged_but_not_serialized() {
        let mut stats = GraphGenerationStats {
            cff_energy_degree_bound_reports: vec![CffEnergyDegreeBoundReport {
                source_kind: CffEnergyBoundSourceKind::PhysicalGraph,
                physical_parent_bounds: vec![(2, 1), (3, 2)],
                assigned_cff_source_bounds: vec![(2, 1), (3, 2)],
            }],
            ..GraphGenerationStats::default()
        };
        stats.merge_in_place(&GraphGenerationStats {
            cff_energy_degree_bound_reports: vec![
                CffEnergyDegreeBoundReport {
                    source_kind: CffEnergyBoundSourceKind::PhysicalGraph,
                    physical_parent_bounds: vec![(2, 1), (3, 2)],
                    assigned_cff_source_bounds: vec![(2, 1), (3, 2)],
                },
                CffEnergyDegreeBoundReport {
                    source_kind: CffEnergyBoundSourceKind::ExactFourD,
                    physical_parent_bounds: vec![(5, 2)],
                    assigned_cff_source_bounds: vec![(9, 1), (10, 1)],
                },
            ],
            ..GraphGenerationStats::default()
        });
        assert_eq!(
            stats.cff_energy_degree_bound_reports,
            vec![
                CffEnergyDegreeBoundReport {
                    source_kind: CffEnergyBoundSourceKind::PhysicalGraph,
                    physical_parent_bounds: vec![(2, 1), (3, 2)],
                    assigned_cff_source_bounds: vec![(2, 1), (3, 2)],
                },
                CffEnergyDegreeBoundReport {
                    source_kind: CffEnergyBoundSourceKind::ExactFourD,
                    physical_parent_bounds: vec![(5, 2)],
                    assigned_cff_source_bounds: vec![(9, 1), (10, 1)],
                },
            ]
        );

        let json = serde_json::to_value(&stats).unwrap();
        assert!(
            json.get("cff_energy_degree_bound_reports").is_none(),
            "transient CFF diagnostics must not change generation_summary.json"
        );
        let decoded: GraphGenerationStats = serde_json::from_value(json).unwrap();
        assert!(decoded.cff_energy_degree_bound_reports.is_empty());
    }
}
