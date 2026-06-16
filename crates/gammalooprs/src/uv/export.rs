use color_eyre::Result;
use eyre::{Context, eyre};
use linnet::half_edge::subgraph::{SubGraphLike, SubSetOps};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    cff::CutCFFIndex,
    graph::{FeynmanGraph, Graph, cuts::CutSet, cuts::ResidueSelector, parse::ToQuoted},
    integrands::process::ProcessIntegrand,
    processes::DotExportSettings,
    settings::global::GenerationSettings,
    uv::{
        UVOrchestrator,
        approx::{CutStructure, OrientationProjection},
        forest::CutForests,
        hedge_poset::Wood as HedgePosetWood,
        wood::CutWoods,
    },
};

pub struct UVForestExportSettings {
    pub computed: bool,
}

pub struct UVForestExport {
    pub graph_name: String,
    pub forest_dot: String,
    pub node_terms: Vec<UVForestNodeTerm>,
}

pub struct UVForestNodeTerm {
    pub forest_index: usize,
    pub node_index: usize,
    pub node_key: String,
    pub term_index: usize,
    pub residue_index: CutCFFIndex,
    pub dot: String,
}

impl UVForestNodeTerm {
    pub fn file_name(&self) -> String {
        let key = sanitize_file_component(&self.node_key);
        let residue = residue_suffix(self.residue_index);
        format!(
            "node_{:03}_{}_term_{:03}_{}.dot",
            self.node_index, key, self.term_index, residue
        )
    }
}

pub(crate) struct UVForestNodeExpression {
    pub forest_index: usize,
    pub node_index: usize,
    pub node_key: String,
    pub term_index: usize,
    pub residue_index: CutCFFIndex,
    pub numerator: Atom,
}

impl ProcessIntegrand {
    pub fn export_uv_forest_graph(
        &self,
        graph_id: usize,
        generation_settings: &GenerationSettings,
        export_settings: &UVForestExportSettings,
    ) -> Result<UVForestExport> {
        match self {
            Self::Amplitude(integrand) => {
                let term = integrand.data.graph_terms.get(graph_id).ok_or_else(|| {
                    eyre!(
                        "Graph id {} is out of range for amplitude integrand {}",
                        graph_id,
                        integrand.data.name
                    )
                })?;
                let cut_structure = CutStructure::empty(&term.graph);
                let post_factor = (Atom::var(crate::utils::GS.pi) * Atom::num(2))
                    .pow(3 * term.graph.get_loop_number() as i64);
                export_graph(
                    &term.graph,
                    cut_structure,
                    term.orientations.iter().cloned().collect(),
                    generation_settings,
                    export_settings,
                    |atom| atom / &post_factor,
                )
            }
            Self::CrossSection(integrand) => {
                let term = integrand.data.graph_terms.get(graph_id).ok_or_else(|| {
                    eyre!(
                        "Graph id {} is out of range for cross-section integrand {}",
                        graph_id,
                        integrand.data.name
                    )
                })?;
                let cuts = term
                    .raised_data
                    .raised_cut_groups
                    .iter()
                    .map(|cuts| CutSet {
                        residue_selector: ResidueSelector {
                            lu_cut: Some(cuts.related_esurface_group.clone()),
                            left_th_cut: None,
                            right_th_cut: None,
                        },
                        union: cuts
                            .cuts
                            .iter()
                            .map(|cut_id| term.cuts[*cut_id].cut.as_subgraph())
                            .reduce(|cut_1, cut_2| cut_1.union(&cut_2))
                            .unwrap_or_else(|| term.graph.empty_subgraph()),
                        canonicalize_external_shifts: false,
                    })
                    .collect();
                let cut_structure = CutStructure { cuts };
                let loop_number = term.graph.cyclotomatic_number(&term.graph.full_filter())
                    - term.graph.initial_state_cut.nedges(&term.graph);
                let loop_3 = loop_number as i64 * 3;
                let lu_prefactor = Atom::var(crate::utils::GS.rescale_star).pow(loop_3)
                    * Atom::var(crate::utils::GS.hfunction_lu_cut)
                    / (Atom::num(2) * Atom::var(crate::utils::GS.pi)).pow(loop_3 - 1);

                export_graph(
                    &term.graph,
                    cut_structure,
                    term.orientations.iter().cloned().collect(),
                    generation_settings,
                    export_settings,
                    |atom| atom * &lu_prefactor,
                )
            }
        }
    }
}

fn export_graph(
    graph: &Graph,
    cut_structure: CutStructure,
    orientations: Vec<
        linnet::half_edge::involution::EdgeVec<linnet::half_edge::involution::Orientation>,
    >,
    generation_settings: &GenerationSettings,
    export_settings: &UVForestExportSettings,
    mut post_process: impl FnMut(Atom) -> Atom,
) -> Result<UVForestExport> {
    if generation_settings.uv.orchestrator == UVOrchestrator::Compare {
        return Err(eyre!(
            "UV forest export does not support uv.orchestrator = compare"
        ));
    }

    let mut forest_dot = String::new();
    let mut node_terms = Vec::new();
    for (forest_index, cut) in cut_structure.cuts.into_iter().enumerate() {
        let single_cut = CutStructure { cuts: vec![cut] };
        let forest_name = format!(
            "uv_{}_forest_{forest_index:03}",
            sanitize_file_component(&graph.name)
        );
        let mut graph = graph.clone();
        match generation_settings.uv.orchestrator {
            UVOrchestrator::LegacyDagForest => export_legacy_forest(
                forest_index,
                &forest_name,
                &mut graph,
                single_cut,
                &orientations,
                generation_settings,
                export_settings,
                &mut post_process,
                &mut forest_dot,
                &mut node_terms,
            )?,
            UVOrchestrator::HedgePoset => export_hedge_poset_forest(
                forest_index,
                &forest_name,
                &mut graph,
                single_cut,
                &orientations,
                generation_settings,
                export_settings,
                &mut post_process,
                &mut forest_dot,
                &mut node_terms,
            )?,
            UVOrchestrator::Compare => unreachable!("compare is rejected before export"),
        }
    }

    Ok(UVForestExport {
        graph_name: graph.name.clone(),
        forest_dot,
        node_terms,
    })
}

#[allow(clippy::too_many_arguments)]
fn export_legacy_forest(
    forest_index: usize,
    forest_name: &str,
    graph: &mut Graph,
    cut_structure: CutStructure,
    orientations: &[linnet::half_edge::involution::EdgeVec<
        linnet::half_edge::involution::Orientation,
    >],
    generation_settings: &GenerationSettings,
    export_settings: &UVForestExportSettings,
    post_process: &mut impl FnMut(Atom) -> Atom,
    forest_dot: &mut String,
    node_terms: &mut Vec<UVForestNodeTerm>,
) -> Result<()> {
    let cut_woods = CutWoods::new(cut_structure, graph, &generation_settings.uv);
    let mut cut_forests = cut_woods.unfold(graph);
    let Some(forest) = cut_forests.forests.first_mut() else {
        return Err(eyre!("Legacy UV exporter produced no forest"));
    };
    forest_dot.push_str(&forest.dot_serialize_for_export(forest_name));
    forest_dot.push('\n');

    if !export_settings.computed {
        return Ok(());
    }

    compute_legacy_forest(graph, &mut cut_forests, orientations, generation_settings)?;
    let forest = cut_forests
        .forests
        .first()
        .expect("legacy forest exists after compute");
    node_terms.extend(
        forest
            .export_node_expressions(forest_index, post_process)?
            .into_iter()
            .map(|term| node_expression_to_dot(graph, forest_name, term))
            .collect::<Result<Vec<_>>>()?,
    );
    Ok(())
}

fn compute_legacy_forest(
    graph: &mut Graph,
    cut_forests: &mut CutForests,
    orientations: &[linnet::half_edge::involution::EdgeVec<
        linnet::half_edge::involution::Orientation,
    >],
    generation_settings: &GenerationSettings,
) -> Result<()> {
    cut_forests.compute(
        graph,
        crate::utils::vakint()?,
        OrientationProjection::new(orientations, &generation_settings.orientation_pattern),
        &generation_settings.uv,
    )
}

#[allow(clippy::too_many_arguments)]
fn export_hedge_poset_forest(
    forest_index: usize,
    forest_name: &str,
    graph: &mut Graph,
    cut_structure: CutStructure,
    orientations: &[linnet::half_edge::involution::EdgeVec<
        linnet::half_edge::involution::Orientation,
    >],
    generation_settings: &GenerationSettings,
    export_settings: &UVForestExportSettings,
    post_process: &mut impl FnMut(Atom) -> Atom,
    forest_dot: &mut String,
    node_terms: &mut Vec<UVForestNodeTerm>,
) -> Result<()> {
    let wood = HedgePosetWood::new(cut_structure, graph, &generation_settings.uv);
    let mut forests = wood.unfold();
    forest_dot.push_str(&name_dot_graph(forests.dot_serialize(), forest_name));
    forest_dot.push('\n');

    if !export_settings.computed {
        return Ok(());
    }

    forests.compute(
        graph,
        crate::utils::vakint()?,
        OrientationProjection::new(orientations, &generation_settings.orientation_pattern),
        &generation_settings.uv,
    )?;
    node_terms.extend(
        forests
            .export_node_expressions(forest_index, post_process)?
            .into_iter()
            .map(|term| node_expression_to_dot(graph, forest_name, term))
            .collect::<Result<Vec<_>>>()?,
    );
    Ok(())
}

fn node_expression_to_dot(
    graph: &Graph,
    forest_name: &str,
    term: UVForestNodeExpression,
) -> Result<UVForestNodeTerm> {
    let graph_name = format!(
        "{}_node_{:03}_term_{:03}",
        forest_name, term.node_index, term.term_index
    );
    let mut dot_graph = graph
        .with_global_numerator_only(graph_name, term.numerator.clone())
        .to_dot_graph_with_settings(&DotExportSettings {
            split_xs_by_initial_states: true,
            output_full_numerator: false,
            ..DotExportSettings::default()
        });
    dot_graph
        .global_data
        .statements
        .insert("forest_name".into(), dot_statement_value(forest_name));
    dot_graph
        .global_data
        .statements
        .insert("forest_index".into(), term.forest_index.to_string());
    dot_graph
        .global_data
        .statements
        .insert("forest_node_index".into(), term.node_index.to_string());
    dot_graph.global_data.statements.insert(
        "forest_node_key".into(),
        dot_statement_value(&term.node_key),
    );
    dot_graph
        .global_data
        .statements
        .insert("forest_term_index".into(), term.term_index.to_string());
    dot_graph.global_data.statements.insert(
        "forest_residue_index".into(),
        dot_statement_value(&residue_suffix(term.residue_index)),
    );
    dot_graph
        .global_data
        .statements
        .insert("full_num".into(), term.numerator.to_quoted());

    let mut dot = Vec::new();
    dot_graph
        .write_io(&mut dot)
        .context("while serializing UV forest node DOT")?;
    Ok(UVForestNodeTerm {
        forest_index: term.forest_index,
        node_index: term.node_index,
        node_key: term.node_key,
        term_index: term.term_index,
        residue_index: term.residue_index,
        dot: String::from_utf8(dot).context("DOT graph serialization was not UTF-8")?,
    })
}

fn name_dot_graph(dot: String, name: &str) -> String {
    dot.replacen("digraph {", &format!("digraph {name} {{"), 1)
        .replacen("digraph Poset {", &format!("digraph {name} {{"), 1)
}

fn escaped_dot_string(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t")
}

pub(crate) fn dot_attr_value(value: &str) -> String {
    format!("\"{}\"", escaped_dot_string(value))
}

fn dot_statement_value(value: &str) -> String {
    escaped_dot_string(value)
}

fn sanitize_file_component(value: &str) -> String {
    let mut sanitized = String::new();
    for c in value.chars() {
        let next = if c.is_ascii_alphanumeric() || c == '-' {
            c
        } else {
            '_'
        };
        if next == '_' {
            if !sanitized.ends_with('_') {
                sanitized.push(next);
            }
        } else {
            sanitized.push(next);
        }
    }
    sanitized.truncate(48);
    sanitized = sanitized.trim_matches('_').to_string();
    if sanitized.is_empty() {
        "key".to_string()
    } else {
        sanitized
    }
}

fn residue_suffix(index: CutCFFIndex) -> String {
    let suffix = index.to_string();
    if suffix.is_empty() {
        "all_none".to_string()
    } else {
        sanitize_file_component(&suffix)
    }
}

#[cfg(test)]
mod tests {
    use super::UVForestNodeTerm;
    use crate::cff::CutCFFIndex;

    #[test]
    fn node_term_file_name_contains_stable_indices_key_and_residue() {
        let term = UVForestNodeTerm {
            forest_index: 2,
            node_index: 4,
            node_key: "{1,2}(3)".to_string(),
            term_index: 7,
            residue_index: CutCFFIndex {
                left_threshold_order: None,
                right_threshold_order: None,
                lu_cut_order: Some(1),
            },
            dot: String::new(),
        };

        assert_eq!(term.file_name(), "node_004_1_2_3_term_007_lu_cut_1.dot");
    }

    #[test]
    fn node_term_file_name_uses_all_none_residue_suffix() {
        let term = UVForestNodeTerm {
            forest_index: 0,
            node_index: 0,
            node_key: "!!!".to_string(),
            term_index: 0,
            residue_index: CutCFFIndex::new_all_none(),
            dot: String::new(),
        };

        assert_eq!(term.file_name(), "node_000_key_term_000_all_none.dot");
    }
}
