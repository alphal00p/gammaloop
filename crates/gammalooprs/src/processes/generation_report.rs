use std::time::Duration;

use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::utils::symbolica_ext::{get_alias_expanded_byte_size, get_all_aliases};

#[derive(Debug, Clone)]
pub struct EvaluatorBuildTimings {
    pub spenso_time: Duration,
    pub symbolica_time: Duration,
    pub atom_byte_size: usize,
    pub atom_alias_expanded_byte_size: usize,
    pub atom_alias_body_byte_size: usize,
    alias_body_atoms: Option<Vec<Atom>>,
}

impl Default for EvaluatorBuildTimings {
    fn default() -> Self {
        Self {
            spenso_time: Duration::default(),
            symbolica_time: Duration::default(),
            atom_byte_size: 0,
            atom_alias_expanded_byte_size: 0,
            atom_alias_body_byte_size: 0,
            alias_body_atoms: Some(Vec::new()),
        }
    }
}

impl EvaluatorBuildTimings {
    pub fn record_alias_atom_size(&mut self, atom: AtomView<'_>) {
        self.atom_byte_size += atom.get_byte_size();
        self.atom_alias_expanded_byte_size += get_alias_expanded_byte_size(atom);
        let aliases = get_all_aliases(atom);
        merge_alias_body_atoms(&mut self.alias_body_atoms, Some(aliases));
        self.atom_alias_body_byte_size = alias_body_byte_size(&self.alias_body_atoms);
    }
}

impl std::ops::AddAssign for EvaluatorBuildTimings {
    fn add_assign(&mut self, other: Self) {
        self.spenso_time += other.spenso_time;
        self.symbolica_time += other.symbolica_time;
        self.atom_byte_size += other.atom_byte_size;
        self.atom_alias_expanded_byte_size += other.atom_alias_expanded_byte_size;
        merge_alias_body_atoms(&mut self.alias_body_atoms, other.alias_body_atoms);
        self.atom_alias_body_byte_size = match &self.alias_body_atoms {
            Some(_) => alias_body_byte_size(&self.alias_body_atoms),
            None => self.atom_alias_body_byte_size + other.atom_alias_body_byte_size,
        };
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphGenerationStats {
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
    #[serde(default)]
    pub evaluator_atom_byte_size: usize,
    #[serde(default)]
    pub evaluator_atom_alias_expanded_byte_size: usize,
    #[serde(default)]
    pub evaluator_atom_alias_body_byte_size: usize,
    #[serde(skip)]
    pub(crate) evaluator_alias_body_atoms: Option<Vec<Atom>>,
}

impl Default for GraphGenerationStats {
    fn default() -> Self {
        Self {
            evaluator_count: 0,
            total_time: Duration::default(),
            evaluator_spenso_time: Duration::default(),
            evaluator_symbolica_time: Duration::default(),
            evaluator_compile_time: Duration::default(),
            evaluator_atom_byte_size: 0,
            evaluator_atom_alias_expanded_byte_size: 0,
            evaluator_atom_alias_body_byte_size: 0,
            evaluator_alias_body_atoms: Some(Vec::new()),
        }
    }
}

impl GraphGenerationStats {
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
        self.evaluator_atom_byte_size += timings.atom_byte_size;
        self.evaluator_atom_alias_expanded_byte_size += timings.atom_alias_expanded_byte_size;
        merge_alias_body_atoms(
            &mut self.evaluator_alias_body_atoms,
            timings.alias_body_atoms,
        );
        self.evaluator_atom_alias_body_byte_size = match &self.evaluator_alias_body_atoms {
            Some(_) => alias_body_byte_size(&self.evaluator_alias_body_atoms),
            None => self.evaluator_atom_alias_body_byte_size + timings.atom_alias_body_byte_size,
        };
    }

    pub fn evaluator_expression_byte_size(&self) -> usize {
        self.evaluator_atom_byte_size + self.evaluator_atom_alias_body_byte_size
    }

    pub fn finalize_alias_size_accounting(&mut self) {
        if let Some(alias_body_atoms) = self.evaluator_alias_body_atoms.take() {
            self.evaluator_atom_alias_body_byte_size =
                alias_body_byte_size_from_atoms(&alias_body_atoms);
        }
    }

    pub fn evaluator_alias_compression_factor(&self) -> Option<f64> {
        let expression_byte_size = self.evaluator_expression_byte_size();
        (expression_byte_size > 0).then_some(
            self.evaluator_atom_alias_expanded_byte_size as f64 / expression_byte_size as f64,
        )
    }

    pub fn merge_in_place(&mut self, other: &Self) {
        self.evaluator_count += other.evaluator_count;
        self.total_time += other.total_time;
        self.evaluator_spenso_time += other.evaluator_spenso_time;
        self.evaluator_symbolica_time += other.evaluator_symbolica_time;
        self.evaluator_compile_time += other.evaluator_compile_time;
        self.evaluator_atom_byte_size += other.evaluator_atom_byte_size;
        self.evaluator_atom_alias_expanded_byte_size +=
            other.evaluator_atom_alias_expanded_byte_size;
        merge_alias_body_atoms(
            &mut self.evaluator_alias_body_atoms,
            other.evaluator_alias_body_atoms.clone(),
        );
        self.evaluator_atom_alias_body_byte_size = match &self.evaluator_alias_body_atoms {
            Some(_) => alias_body_byte_size(&self.evaluator_alias_body_atoms),
            None => {
                self.evaluator_atom_alias_body_byte_size + other.evaluator_atom_alias_body_byte_size
            }
        };
    }
}

fn merge_alias_body_atoms(dst: &mut Option<Vec<Atom>>, src: Option<Vec<Atom>>) {
    match (dst, src) {
        (Some(dst), Some(src)) => {
            for atom in src {
                if !dst
                    .iter()
                    .any(|existing| existing.as_atom_view() == atom.as_atom_view())
                {
                    dst.push(atom);
                }
            }
        }
        (dst @ Some(_), None) => *dst = None,
        (None, Some(_)) => {}
        (None, None) => {}
    }
}

fn alias_body_byte_size(alias_body_atoms: &Option<Vec<Atom>>) -> usize {
    alias_body_atoms
        .as_ref()
        .map(|atoms| alias_body_byte_size_from_atoms(atoms))
        .unwrap_or(0)
}

fn alias_body_byte_size_from_atoms(alias_body_atoms: &[Atom]) -> usize {
    alias_body_atoms
        .iter()
        .map(|atom| atom.as_atom_view().get_byte_size())
        .sum()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NamedGraphGenerationReport {
    pub integrand_name: String,
    pub graph_name: String,
    pub stats: GraphGenerationStats,
}

impl NamedGraphGenerationReport {
    pub fn finalize_alias_size_accounting(&mut self) {
        self.stats.finalize_alias_size_accounting();
    }
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

    pub fn finalize_alias_size_accounting(&mut self) {
        self.stats.finalize_alias_size_accounting();
    }
}

pub fn finalize_generated_graph_reports_alias_size_accounting(
    reports: &mut [GeneratedGraphReport],
) {
    for report in reports {
        report.finalize_alias_size_accounting();
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
