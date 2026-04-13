use std::time::Duration;

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GraphGenerationStats {
    pub evaluator_count: usize,
    pub total_time: Duration,
    pub evaluator_build_time: Duration,
    pub evaluator_compile_time: Duration,
}

impl GraphGenerationStats {
    pub fn expression_build_time(&self) -> Duration {
        self.total_time
            .saturating_sub(self.evaluator_build_time)
            .saturating_sub(self.evaluator_compile_time)
    }

    pub fn merge_in_place(&mut self, other: &Self) {
        self.evaluator_count += other.evaluator_count;
        self.total_time += other.total_time;
        self.evaluator_build_time += other.evaluator_build_time;
        self.evaluator_compile_time += other.evaluator_compile_time;
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
