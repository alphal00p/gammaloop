use std::sync::Arc;

use color_eyre::Result;
use eyre::{Context, eyre};
use linnet::half_edge::subgraph::SubSetOps;
use linnet::parser::DotGraph;
use serde::{Deserialize, Serialize};
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use symbolica::atom::{Atom, AtomCore};

use crate::{
    cff::CutCFFIndex,
    graph::{
        Graph,
        cuts::{CutSet, ResidueSelector},
    },
    integrands::process::ProcessIntegrand,
    processes::{DotExportSettings, ProcessCollection, process::Process},
    settings::global::GenerationSettings,
    uv::{
        UVOrchestrator,
        approx::{
            CutStructure, OrientationProjection,
            local_3d::Local3DProjectionPath,
            local_4d::{Local4dBranchProvenance, Local4dComponentProvenance},
        },
        forest::CutForests,
        hedge_poset::Wood as HedgePosetWood,
        orchestrator::validate_local_counterterm_schemes,
        settings::FinalIntegrandDimension,
        wood::CutWoods,
    },
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct UVForestNodeProvenance {
    pub parent_keys: Vec<String>,
    pub local: UVForestLocalProvenance,
}

/// The provenance representation must match the numerator exported beside it.
/// A direct 3D export therefore records its CFF projection paths rather than
/// forcing construction of an unrelated 4D counterterm atom.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "representation")]
pub(crate) enum UVForestLocalProvenance {
    #[serde(rename = "4d")]
    FourD {
        components: Vec<Local4dComponentProvenance>,
        branches: Vec<Local4dBranchProvenance>,
    },
    #[serde(rename = "3d")]
    ThreeD {
        projection_paths: Vec<Local3DProjectionPath>,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct UVForestProvenance {
    pub node: UVForestNodeProvenance,
}

pub struct UVForestExportSettings {
    /// Whether to evaluate forest nodes. Computed exports consume the
    /// orientations already selected while generating the graph term;
    /// structure-only exports do not evaluate orientation-dependent atoms.
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

    /// Return the atom-free provenance embedded in this computed Hedge DOT.
    ///
    /// Parsing remains encapsulated here so callers do not need a direct DOT
    /// parser dependency. Legacy and structure-only exports return `None`.
    pub fn forest_provenance(&self) -> Result<Option<serde_json::Value>> {
        let dot: DotGraph = DotGraph::from_string(&self.dot)
            .map_err(|error| eyre!("failed to parse UV forest node DOT: {error}"))?;
        let Some(json) = dot.global_data.statements.get("forest_provenance") else {
            return Ok(None);
        };
        Ok(Some(
            serde_json::from_str(json).context("while parsing UV forest provenance JSON")?,
        ))
    }
}

pub(crate) struct UVForestNodeExpression {
    pub forest_index: usize,
    pub node_index: usize,
    pub node_key: String,
    pub term_index: usize,
    pub residue_index: CutCFFIndex,
    pub numerator: Atom,
    pub forest_provenance: Option<Arc<UVForestProvenance>>,
}

impl Process {
    pub fn export_uv_forest_graph(
        &self,
        integrand_name: &str,
        graph_id: usize,
        settings: &UVForestExportSettings,
    ) -> Result<UVForestExport> {
        let generation_settings = &self
            .settings_history
            .as_ref()
            .ok_or_else(|| {
                eyre!(
                    "Cannot export UV forests for process {} without generation settings history",
                    self.definition.folder_name
                )
            })?
            .generation;
        let resolved = self.get_integrand(integrand_name)?;
        let integrand = resolved.require_generated()?;
        let source = if settings.computed {
            let (graph, expression) = match &self.collection {
                ProcessCollection::Amplitudes(amplitudes) => {
                    let source = amplitudes[&resolved.canonical_name]
                        .graphs
                        .get(graph_id)
                        .ok_or_else(|| eyre!("Missing source amplitude graph {graph_id}"))?;
                    (&source.graph, source.derived_data.cff_expression.as_ref())
                }
                ProcessCollection::CrossSections(cross_sections) => {
                    let source = cross_sections[&resolved.canonical_name]
                        .supergraphs
                        .get(graph_id)
                        .ok_or_else(|| eyre!("Missing source cross-section graph {graph_id}"))?;
                    (
                        &source.graph,
                        source.derived_data.global_cff_expression.as_ref(),
                    )
                }
            };
            // Generation and persistent selection preserve graph order;
            // reject stale runtime metadata instead of pairing a different
            // graph with this stored production residue map.
            if integrand.graph_name_by_id(graph_id) != Some(graph.name.as_str()) {
                return Err(eyre!(
                    "Source/runtime graph mismatch for computed UV forest export at id {graph_id}"
                ));
            }
            Some((
                graph,
                expression.ok_or_else(|| {
                    eyre!(
                        "Graph {} has no stored production CFF for computed UV forest export",
                        graph.name
                    )
                })?,
            ))
        } else {
            None
        };
        let cff_options = source
            .map(|(graph, _)| graph.production_cff_3d_expression_options(generation_settings))
            .transpose()?;
        let orientation = source
            .zip(cff_options.as_ref())
            .map(|((_, expression), options)| {
                OrientationProjection::exact_expression(
                    expression,
                    options,
                    &generation_settings.orientation_pattern,
                    generation_settings.explicit_orientation_sum_only,
                )
            });
        match integrand {
            ProcessIntegrand::Amplitude(integrand) => {
                let term = integrand.data.graph_terms.get(graph_id).ok_or_else(|| {
                    eyre!(
                        "Graph id {} is out of range for amplitude integrand {}",
                        graph_id,
                        integrand.data.name
                    )
                })?;
                let cut_structure = CutStructure::empty(&term.graph);
                export_graph(
                    &term.graph,
                    cut_structure,
                    orientation,
                    generation_settings,
                    settings,
                    |_, atom| atom,
                )
            }
            ProcessIntegrand::CrossSection(integrand) => {
                let term = integrand.data.graph_terms.get(graph_id).ok_or_else(|| {
                    eyre!(
                        "Graph id {} is out of range for cross-section integrand {}",
                        graph_id,
                        integrand.data.name
                    )
                })?;
                let cuts = term
                    .cut_group_data
                    .cut_groups
                    .iter()
                    .map(|cuts| CutSet {
                        residue_selector: ResidueSelector {
                            lu: Some(cuts.lu_cut_selection(&term.graph, &term.cuts)),
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
                // Finalized UV nodes already include the CFF spatial measure.
                // Reuse the physical cut-group conversion, as production does.
                let lu_prefactors = term
                    .cut_group_data
                    .cut_groups
                    .iter()
                    .map(|group| group.lu_prefactor(&term.graph, &term.cuts))
                    .collect::<Result<Vec<_>>>()?;

                export_graph(
                    &term.graph,
                    cut_structure,
                    orientation,
                    generation_settings,
                    settings,
                    |forest_index, atom| atom * &lu_prefactors[forest_index],
                )
            }
        }
    }
}

fn export_graph(
    graph: &Graph,
    cut_structure: CutStructure,
    orientation: Option<OrientationProjection<'_>>,
    generation_settings: &GenerationSettings,
    export_settings: &UVForestExportSettings,
    mut post_process: impl FnMut(usize, Atom) -> Atom,
) -> Result<UVForestExport> {
    if generation_settings.uv.orchestrator == UVOrchestrator::Compare {
        return Err(eyre!(
            "UV forest export does not support uv.orchestrator = compare"
        ));
    }
    validate_local_counterterm_schemes(graph, &cut_structure, &generation_settings.uv)?;

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
                orientation,
                generation_settings,
                export_settings,
                &mut |atom| post_process(forest_index, atom),
                &mut forest_dot,
                &mut node_terms,
            )?,
            UVOrchestrator::HedgePoset => export_hedge_poset_forest(
                forest_index,
                &forest_name,
                &mut graph,
                single_cut,
                orientation,
                generation_settings,
                export_settings,
                &mut |atom| post_process(forest_index, atom),
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
    orientation: Option<OrientationProjection<'_>>,
    generation_settings: &GenerationSettings,
    export_settings: &UVForestExportSettings,
    post_process: &mut impl FnMut(Atom) -> Atom,
    forest_dot: &mut String,
    node_terms: &mut Vec<UVForestNodeTerm>,
) -> Result<()> {
    let cut_woods = CutWoods::new(cut_structure, graph, &generation_settings.uv)?;
    let mut cut_forests = cut_woods.unfold(graph);
    let Some(forest) = cut_forests.forests.first_mut() else {
        return Err(eyre!("Legacy UV exporter produced no forest"));
    };
    forest_dot.push_str(&forest.dot_serialize_for_export(forest_name));
    forest_dot.push('\n');

    if !export_settings.computed {
        return Ok(());
    }

    if matches!(
        generation_settings.uv.final_integrand,
        FinalIntegrandDimension::FourD
    ) {
        return Err(eyre!(
            "computed 4D UV forest export is not supported by the legacy DAG forest for graph '{}' (forest {forest_index}); use the HedgePoset orchestrator",
            graph.name,
        ));
    }
    let orientation = orientation.ok_or_else(|| {
        eyre!("Computed UV forest export requires its stored production CFF expression")
    })?;
    compute_legacy_forest(graph, &mut cut_forests, orientation, generation_settings)?;
    let forest = cut_forests
        .forests
        .first()
        .expect("legacy forest exists after compute");
    node_terms.extend(
        forest
            .export_node_expressions(graph, forest_index, post_process)?
            .into_iter()
            .map(|term| node_expression_to_dot(graph, forest_name, term))
            .collect::<Result<Vec<_>>>()?,
    );
    Ok(())
}

fn compute_legacy_forest(
    graph: &mut Graph,
    cut_forests: &mut CutForests,
    orientation: OrientationProjection<'_>,
    generation_settings: &GenerationSettings,
) -> Result<()> {
    cut_forests.compute(
        graph,
        crate::utils::vakint()?,
        orientation,
        &generation_settings.uv,
    )
}

#[allow(clippy::too_many_arguments)]
fn export_hedge_poset_forest(
    forest_index: usize,
    forest_name: &str,
    graph: &mut Graph,
    cut_structure: CutStructure,
    orientation: Option<OrientationProjection<'_>>,
    generation_settings: &GenerationSettings,
    export_settings: &UVForestExportSettings,
    post_process: &mut impl FnMut(Atom) -> Atom,
    forest_dot: &mut String,
    node_terms: &mut Vec<UVForestNodeTerm>,
) -> Result<()> {
    let wood = HedgePosetWood::new(cut_structure, graph, &generation_settings.uv)?;
    let mut forests = wood.unfold();
    forest_dot.push_str(&name_dot_graph(forests.dot_serialize(), forest_name));
    forest_dot.push('\n');

    if !export_settings.computed {
        return Ok(());
    }

    let expressions = match &generation_settings.uv.final_integrand {
        FinalIntegrandDimension::FourD => {
            forests.integrate(graph, crate::utils::vakint()?, &generation_settings.uv)?;
            forests.export_4d_node_expressions(forest_index)?
        }
        FinalIntegrandDimension::ThreeD => {
            // A computed 3D export carries provenance from the direct
            // orientation-local CFF projections themselves. Do not construct
            // unrelated 4D atoms merely to annotate this diagnostic output.
            forests.compute(
                graph,
                crate::utils::vakint()?,
                orientation.ok_or_else(|| {
                    eyre!("Computed UV forest export requires its stored production CFF expression")
                })?,
                &generation_settings.uv,
            )?;
            forests.export_node_expressions(forest_index, post_process)?
        }
    };
    node_terms.extend(
        expressions
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
        .insert("forest_name".into(), forest_name.to_owned());
    dot_graph
        .global_data
        .statements
        .insert("forest_index".into(), term.forest_index.to_string());
    dot_graph
        .global_data
        .statements
        .insert("forest_node_index".into(), term.node_index.to_string());
    dot_graph
        .global_data
        .statements
        .insert("forest_node_key".into(), term.node_key.clone());
    dot_graph
        .global_data
        .statements
        .insert("forest_term_index".into(), term.term_index.to_string());
    dot_graph.global_data.statements.insert(
        "forest_residue_index".into(),
        residue_suffix(term.residue_index),
    );
    if let Some(provenance) = term.forest_provenance.as_deref() {
        let provenance = serde_json::to_string(provenance)
            .context("while serializing UV forest provenance as JSON")?;
        dot_graph
            .global_data
            .statements
            .insert("forest_provenance".into(), provenance);
    }
    let full_num = term
        .numerator
        .printer(SpensoPrintSettings::typst_options())
        .to_string();
    dot_graph
        .global_data
        .statements
        .insert("full_num".into(), full_num);

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

pub(crate) fn sanitize_file_component(value: &str) -> String {
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
    use std::sync::{Arc, Once, OnceLock};

    use linnet::parser::DotGraph;
    use symbolica::{
        atom::{Atom, AtomCore},
        function, symbol,
    };

    use super::{
        UVForestExportSettings, UVForestLocalProvenance, UVForestNodeExpression,
        UVForestNodeProvenance, UVForestNodeTerm, UVForestProvenance, export_graph,
        node_expression_to_dot, sanitize_file_component,
    };
    use crate::{
        cff::CutCFFIndex,
        dot,
        graph::{FeynmanGraph, Graph, parse::IntoGraph},
        initialisation::test_initialise,
        settings::global::GenerationSettings,
        utils::GS,
        uv::{
            ApproximationType, UVOrchestrator, UVgenerationSettings,
            approx::{
                CutStructure, OrientationProjection,
                local_4d::{
                    Local4dBranch, Local4dBranchProvenance, Local4dComponentProvenance,
                    Local4dRouteSignature,
                },
            },
            settings::FinalIntegrandDimension,
        },
    };

    static EXPORT_TEST_INIT: Once = Once::new();
    static SM_MODEL: OnceLock<crate::model::Model> = OnceLock::new();
    static SCALAR_MODEL: OnceLock<crate::model::Model> = OnceLock::new();

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

    #[test]
    fn file_components_cannot_escape_the_export_directory() {
        assert_eq!(sanitize_file_component("../result"), "result");
        assert_eq!(sanitize_file_component("/tmp/result"), "tmp_result");
        assert_eq!(sanitize_file_component("a/b"), "a_b");
        assert_eq!(sanitize_file_component(""), "key");
    }

    #[test]
    fn computed_node_full_numerator_is_a_spenso_aware_typst_fragment() {
        EXPORT_TEST_INIT.call_once(|| test_initialise().unwrap());
        let graph: Graph = dot!(digraph G {
                ext [style=invis]
                node [num=1]
                ext -> A
                C -> A
                A -> D
                D -> B
                B -> C
                C -> D
                B -> ext
            }, SM_MODEL.get_or_init(|| crate::utils::load_generic_model("sm")))
        .unwrap();
        let term = UVForestNodeExpression {
            forest_index: 0,
            node_index: 1,
            node_key: "{1}".to_string(),
            term_index: 2,
            residue_index: CutCFFIndex::new_all_none(),
            numerator: function!(GS.uv_truncate, Atom::num(1)),
            forest_provenance: None,
        };

        let exported = node_expression_to_dot(&graph, "uv_typst", term).unwrap();

        assert!(exported.dot.contains(r#"full_num = "op(\"Tr\")(1)";"#));
        assert!(!exported.dot.contains("gammalooprs::Truncate"));
        assert!(!exported.dot.contains("forest_provenance"));
    }

    #[test]
    fn computed_node_full_numerator_distinguishes_physical_and_uv_mass_denominators() {
        EXPORT_TEST_INIT.call_once(|| test_initialise().unwrap());
        let graph: Graph = dot!(digraph G {
                ext [style=invis]
                node [num=1]
                ext -> A
                A -> A
                A -> ext
            }, SM_MODEL.get_or_init(|| crate::utils::load_generic_model("sm")))
        .unwrap();
        let momentum_squared = Atom::var(symbol!("inspect_momentum_squared"));
        let physical_mass = Atom::var(symbol!("inspect_physical_mass"));
        let uv_expansion_mass = Atom::var(GS.m_uv_expansion);
        let uv_vacuum_mass = Atom::var(GS.m_uv_vacuum);
        let physical_denominator = function!(
            GS.den,
            1,
            momentum_squared.clone(),
            physical_mass.clone().pow(2),
            &momentum_squared - physical_mass.clone().pow(2)
        );
        let terminal_uv_denominator = function!(
            GS.den,
            2,
            momentum_squared.clone(),
            uv_expansion_mass.clone().pow(2),
            momentum_squared - uv_vacuum_mass.clone().pow(2)
        );
        let term = UVForestNodeExpression {
            forest_index: 0,
            node_index: 1,
            node_key: "{1}".to_string(),
            term_index: 0,
            residue_index: CutCFFIndex::new_all_none(),
            numerator: physical_denominator * terminal_uv_denominator,
            forest_provenance: None,
        };

        let exported = node_expression_to_dot(&graph, "uv_mass_provenance", term).unwrap();
        let parsed: DotGraph = DotGraph::from_string(&exported.dot).unwrap();
        let full_num = parsed
            .global_data
            .statements
            .get("full_num")
            .expect("computed inspect DOT must contain its full numerator");

        assert_eq!(
            full_num.matches("denom").count(),
            2,
            "inspect output must retain the physical and terminal-UV denominators separately: {full_num}"
        );
        assert!(
            full_num.contains("inspect_physical_mass"),
            "inspect output lost the physical-mass denominator: {full_num}"
        );
        assert!(
            full_num.contains("UVexp") && full_num.contains("UV"),
            "inspect output lost the distinct UV expansion and vacuum-mass roles: {full_num}"
        );
    }

    #[test]
    fn computed_hedge_4d_export_has_one_all_none_term_per_node() {
        EXPORT_TEST_INIT.call_once(|| test_initialise().unwrap());
        let graph: Graph = include_str!("../../../../tests/resources/graphs/scalar_bubble.dot")
            .into_graph(SCALAR_MODEL.get_or_init(|| crate::utils::load_generic_model("scalars")))
            .unwrap();
        let generation_settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                final_integrand: FinalIntegrandDimension::FourD,
                orchestrator: UVOrchestrator::HedgePoset,
                ..Default::default()
            },
            ..Default::default()
        };
        let export = export_graph(
            &graph,
            CutStructure::empty(&graph),
            None,
            &generation_settings,
            &UVForestExportSettings { computed: true },
            |_, _| panic!("a computed 4D export must not apply 3D post-processing"),
        )
        .unwrap();

        assert!(!export.node_terms.is_empty());
        let mut nodes = std::collections::BTreeSet::new();
        for term in &export.node_terms {
            assert!(nodes.insert((term.forest_index, term.node_index, &term.node_key)));
            assert_eq!(term.term_index, 0);
            assert_eq!(term.residue_index, CutCFFIndex::new_all_none());
            let provenance = term
                .forest_provenance()
                .unwrap()
                .expect("computed Hedge export must retain 4D provenance");
            assert!(provenance.get("local_3d").is_none());
            assert_eq!(provenance.as_object().unwrap().len(), 1);
            assert!(provenance["node"].get("components").is_none());
            assert!(provenance["node"].get("local_4d_branches").is_none());
            assert_eq!(provenance["node"]["local"]["representation"], "4d");
            assert!(provenance["node"]["local"]["components"].is_array());
            assert!(provenance["node"]["local"]["branches"].is_array());
        }
    }

    #[test]
    fn computed_hedge_3d_export_retains_direct_projection_provenance() {
        EXPORT_TEST_INIT.call_once(|| test_initialise().unwrap());
        let mut graph: Graph = include_str!("../../../../tests/resources/graphs/scalar_bubble.dot")
            .into_graph(SCALAR_MODEL.get_or_init(|| crate::utils::load_generic_model("scalars")))
            .unwrap();
        let generation_settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                final_integrand: FinalIntegrandDimension::ThreeD,
                orchestrator: UVOrchestrator::HedgePoset,
                ..Default::default()
            },
            ..Default::default()
        };
        let options = graph
            .production_cff_3d_expression_options(&generation_settings)
            .unwrap();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        let production = graph
            .generate_3d_expression_for_integrand(&[], &canonization, &options, Some(&numerator))
            .unwrap();
        let orientation = OrientationProjection::exact_expression(
            &production,
            &options,
            &generation_settings.orientation_pattern,
            generation_settings.explicit_orientation_sum_only,
        );
        let export = export_graph(
            &graph,
            CutStructure::empty(&graph),
            Some(orientation),
            &generation_settings,
            &UVForestExportSettings { computed: true },
            |_, atom| atom,
        )
        .unwrap();

        assert!(!export.node_terms.is_empty());
        let provenances = export
            .node_terms
            .iter()
            .map(|term| {
                term.forest_provenance()
                    .unwrap()
                    .expect("computed 3D Hedge export must retain direct projection provenance")
            })
            .collect::<Vec<_>>();
        assert!(provenances.iter().all(|provenance| {
            provenance["node"].get("components").is_none()
                && provenance["node"].get("local_4d_branches").is_none()
                && provenance["node"]["local"]["representation"] == "3d"
                && provenance["node"]["local"]["projection_paths"].is_array()
        }));
        assert!(
            provenances
                .iter()
                .any(|provenance| provenance["node"]["local"]["projection_paths"]
                    .as_array()
                    .is_some_and(|paths| paths.iter().any(|path| path["steps"]
                        .as_array()
                        .is_some_and(|steps| !steps.is_empty())))),
            "a non-root computed 3D term must report its actual CFF projection path"
        );
        for step in provenances
            .iter()
            .flat_map(|provenance| {
                provenance["node"]["local"]["projection_paths"]
                    .as_array()
                    .expect("projection paths were validated above")
            })
            .flat_map(|path| {
                path["steps"]
                    .as_array()
                    .expect("a direct projection path must contain a steps array")
            })
        {
            assert!(step["current_component"].is_string());
            assert!(step["given_component"].is_string());
            assert!(step["active_subgraph"].is_string());
            assert!(step["rescaled_subgraph"].is_string());
            assert!(step["canonical_route_loop_edges"].is_array());
            assert!(step["route_loop_edges"].is_array());
            assert_eq!(step["conceptual_branches"][0]["branch"], "U");
            assert_eq!(
                step["conceptual_branches"][0]["coefficient"].as_i64(),
                Some(1)
            );
            assert_eq!(step["materialization"], "direct_u");
        }
    }

    #[test]
    fn computed_legacy_4d_export_is_rejected_contextually() {
        EXPORT_TEST_INIT.call_once(|| test_initialise().unwrap());
        let graph: Graph = include_str!("../../../../tests/resources/graphs/scalar_bubble.dot")
            .into_graph(SCALAR_MODEL.get_or_init(|| crate::utils::load_generic_model("scalars")))
            .unwrap();
        let generation_settings = GenerationSettings {
            uv: UVgenerationSettings {
                final_integrand: FinalIntegrandDimension::FourD,
                orchestrator: UVOrchestrator::LegacyDagForest,
                ..Default::default()
            },
            ..Default::default()
        };
        let error = export_graph(
            &graph,
            CutStructure::empty(&graph),
            None,
            &generation_settings,
            &UVForestExportSettings { computed: true },
            |_, atom| atom,
        )
        .err()
        .expect("legacy computed 4D export must fail before expression generation");

        let message = error.to_string();
        assert!(message.contains("computed 4D UV forest export"));
        assert!(message.contains("legacy DAG forest"));
        assert!(message.contains("graph 'bubble'"));
        assert!(message.contains("forest 0"));
        assert!(message.contains("HedgePoset"));
    }

    #[test]
    fn computed_hedge_provenance_survives_json_and_dot_round_trip() {
        EXPORT_TEST_INIT.call_once(|| test_initialise().unwrap());
        let graph: Graph = dot!(digraph G {
                ext [style=invis]
                node [num=1]
                ext -> A
                A -> A
                A -> ext
            }, SM_MODEL.get_or_init(|| crate::utils::load_generic_model("sm")))
        .unwrap();
        let provenance = Arc::new(UVForestProvenance {
            node: UVForestNodeProvenance {
                parent_keys: vec!["parent:{1}".to_string()],
                local: UVForestLocalProvenance::FourD {
                    components: vec![Local4dComponentProvenance {
                        component: "12".to_string(),
                        scheme: ApproximationType::IR,
                        dod: 2,
                        route_loop_edges: vec![1],
                        route_external_edges: vec![0, 2],
                        route_signatures: vec![Local4dRouteSignature {
                            edge: 1,
                            loop_signature: "+".to_string(),
                            external_signature: "00".to_string(),
                        }],
                    }],
                    branches: vec![
                        Local4dBranchProvenance {
                            branch: Local4dBranch::U,
                            byte_size: 11,
                        },
                        Local4dBranchProvenance {
                            branch: Local4dBranch::S,
                            byte_size: 7,
                        },
                        Local4dBranchProvenance {
                            branch: Local4dBranch::US,
                            byte_size: 5,
                        },
                        Local4dBranchProvenance {
                            branch: Local4dBranch::Combined,
                            byte_size: 13,
                        },
                    ],
                },
            },
        });
        let term = UVForestNodeExpression {
            forest_index: 3,
            node_index: 4,
            node_key: "{1};{2}".to_string(),
            term_index: 5,
            residue_index: CutCFFIndex::new_all_none(),
            numerator: Atom::one(),
            forest_provenance: Some(provenance.clone()),
        };
        let exported = node_expression_to_dot(&graph, "uv_provenance", term).unwrap();
        assert_eq!(
            exported.forest_provenance().unwrap(),
            Some(serde_json::to_value(provenance.as_ref()).unwrap())
        );
        let reparsed: DotGraph = DotGraph::from_string(&exported.dot).unwrap();
        let json = &reparsed.global_data.statements["forest_provenance"];
        let round_trip: UVForestProvenance = serde_json::from_str(json).unwrap();

        assert_eq!(round_trip, *provenance);
        assert_eq!(reparsed.global_data.statements["forest_index"], "3");
        assert_eq!(reparsed.global_data.statements["forest_node_index"], "4");
        assert_eq!(
            reparsed.global_data.statements["forest_node_key"],
            "{1};{2}"
        );
        assert_eq!(reparsed.global_data.statements["forest_term_index"], "5");
        assert_eq!(
            reparsed.global_data.statements["forest_residue_index"],
            "all_none"
        );
    }

    #[test]
    fn computed_exports_match_physical_production_normalization() -> color_eyre::Result<()> {
        use crate::{
            feyngen::GenerationType,
            processes::{Process, ProcessCollection, ProcessDefinition},
            settings::{GlobalSettings, RuntimeSettings},
            utils::{W_, load_generic_model},
            uv::marker::UvMarker,
        };
        use linnet::half_edge::{involution::EdgeIndex, subgraph::SubSetOps};
        use symbolica::atom::{AtomCore, AtomView};

        test_initialise()?;
        let model = load_generic_model("scalars");
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
        let runtime = RuntimeSettings::default();
        for orchestrator in [UVOrchestrator::LegacyDagForest, UVOrchestrator::HedgePoset] {
            // One amplitude contour and two Born cut multiplicities expose both
            // a lost cut i and any second insertion of the spatial measure.
            for multiplicity in [0, 2, 3] {
                let (kind, source) = if multiplicity == 0 {
                    (
                        GenerationType::Amplitude,
                        include_str!(concat!(
                            env!("CARGO_MANIFEST_DIR"),
                            "/../../tests/resources/graphs/scalar_bubble.dot"
                        ))
                        .to_string(),
                    )
                } else {
                    let mut source = String::from(
                        "digraph scalar_born { ext [style=invis]; \
                         ext -> left [particle=scalar_2,is_cut=0]; \
                         right -> ext [particle=scalar_2,is_cut=0];",
                    );
                    for _ in 0..multiplicity {
                        source.push_str("left -> right [particle=scalar_0];");
                    }
                    source.push('}');
                    (GenerationType::CrossSection, source)
                };
                let graphs = Graph::from_string(&source, &model)?;
                let definition = ProcessDefinition::from_graph_list(&graphs, kind, &model)?;
                let mut process = Process::from_graph_list(
                    "export_phase".into(),
                    "default".into(),
                    graphs,
                    kind,
                    Some(definition),
                    None,
                    &model,
                )?;
                let mut settings = GlobalSettings::default();
                settings.generation.uv.orchestrator = orchestrator;
                settings.generation.uv.subtract_uv = false;
                settings.generation.uv.generate_integrated = false;
                settings.generation.threshold_subtraction.enable_thresholds = false;
                settings.generation.explicit_orientation_sum_only = true;
                settings.generation.evaluator.compile = false;
                process.preprocess(&model, &settings, &(&runtime).into(), &pool)?;
                process.generate_integrands(&model, &settings, (&runtime).into(), &pool)?;
                let (graph, expected) = match &process.collection {
                    ProcessCollection::Amplitudes(amplitudes) => {
                        let graph = &amplitudes["default"].graphs[0];
                        (&graph.graph, graph.derived_data.resolved_integrand()?)
                    }
                    ProcessCollection::CrossSections(cross_sections) => {
                        let graph = &cross_sections["default"].supergraphs[0];
                        (
                            &graph.graph,
                            graph
                                .derived_data
                                .cut_paramatric_integrand
                                .iter()
                                .map(|integrand| integrand.integrands.resolved())
                                .collect::<color_eyre::Result<Vec<_>>>()?
                                .iter()
                                .flat_map(|integrands| integrands.iter())
                                .map(|(_, atom)| atom.clone())
                                .sum::<Atom>(),
                        )
                    }
                };
                let exported = process.export_uv_forest_graph(
                    "default",
                    0,
                    &UVForestExportSettings { computed: true },
                )?;
                let split = graph
                    .iter_edges_of(
                        &graph
                            .full_filter()
                            .subtract(&graph.initial_state_cut)
                            .subtract(&graph.tree_edges),
                    )
                    .map(|(_, edge, _)| GS.split_mom_pattern_simple(edge))
                    .collect::<Vec<_>>();
                let mut actual = Atom::Zero;
                for term in &exported.node_terms {
                    for parsed in Graph::from_string(&term.dot, &model)? {
                        actual += UvMarker::new(&settings.generation.uv).finish(
                            &parsed
                                .global_prefactor
                                .num
                                .replace_multiple(&split)
                                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                                .with(W_.d_),
                        );
                    }
                }
                // Compare the stored trees first. If their factorizations differ,
                // evaluate those trees directly at exact nonsingular points.
                // These signed checks certify the sampled values, not a symbolic
                // identity; neither numerator is expanded or polynomialized.
                if actual != expected {
                    let mut inputs = actual
                        .get_all_indeterminates(false)
                        .into_iter()
                        .chain(expected.get_all_indeterminates(false))
                        .map(|input| input.to_owned())
                        .collect::<Vec<_>>();
                    inputs.sort();
                    inputs.dedup();
                    for base in [3, 5, 7] {
                        let mut actual_value = actual.clone();
                        let mut expected_value = expected.clone();
                        for (index, input) in inputs.iter().enumerate() {
                            let value = Atom::num(base).pow(index as i64 + 1);
                            actual_value = actual_value.replace(input.clone()).with(value.clone());
                            expected_value = expected_value.replace(input.clone()).with(value);
                        }
                        assert!(matches!(actual_value.as_view(), AtomView::Num(_)));
                        assert!(matches!(expected_value.as_view(), AtomView::Num(_)));
                        assert!(!expected_value.is_zero());
                        assert_eq!(
                            actual_value, expected_value,
                            "{kind}/{orchestrator}: computed export differs from production at base {base}"
                        );
                    }
                }
                if multiplicity != 0 {
                    let [coupling, expression] = model.get_coupling("SCALAR_COUPLING").rep_rule();
                    let lambda: Atom = model.get_parameter("lam").name.into();
                    // The runtime tree-denominator function is the identity;
                    // these contact Born graphs have an empty product inside it.
                    let mut physical = actual
                        .replace(GS.wrap_tree_denoms(Atom::one()))
                        .with(Atom::one())
                        .replace(coupling)
                        .with(expression)
                        .replace(lambda)
                        .with(Atom::one())
                        .replace(GS.rescale_star)
                        .with(Atom::one())
                        .replace(GS.hfunction_lu_cut)
                        .with(Atom::one());
                    for edge in 0..graph.n_edges() {
                        physical = physical.replace(GS.ose(EdgeIndex(edge))).with(Atom::one());
                    }
                    // With all cut energies set to one, |A|²=1 and no graph
                    // symmetry factor, dPhi_n has residue coefficient
                    // 2pi / [(2pi)^(3(n-1)) * 2^n]. The LU derivative is separate.
                    let two_pi = Atom::num(2) * Atom::var(GS.pi);
                    let born = &two_pi
                        / (two_pi.pow(3 * (multiplicity - 1)) * Atom::num(2).pow(multiplicity));
                    assert_eq!(
                        physical, born,
                        "{orchestrator}: scalar Born export has the wrong phase or normalization"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn computed_export_retains_nonempty_forests_in_all_routes_and_orchestrators()
    -> color_eyre::Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph computed_direct_forest {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let mut generation = GenerationSettings::default();
        assert!(generation.uv.subtract_uv);
        assert!(generation.uv.generate_integrated);
        assert_eq!(
            generation.uv.final_integrand,
            crate::uv::settings::FinalIntegrandDimension::ThreeD
        );
        let options = graph.production_cff_3d_expression_options(&generation)?;
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&Atom::one()),
        )?;
        assert!(!production.expression.orientations.is_empty());

        for orchestrator in [UVOrchestrator::LegacyDagForest, UVOrchestrator::HedgePoset] {
            generation.uv.orchestrator = orchestrator;
            for (explicit_sum, projected) in [(false, false), (true, false), (true, true)] {
                generation.explicit_orientation_sum_only = explicit_sum;
                generation.uv.local_uv_cts_from_expanded_4d_integrands = projected;
                let projection = OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &generation.orientation_pattern,
                    generation.explicit_orientation_sum_only,
                );
                let exported = export_graph(
                    &graph,
                    CutStructure::empty(&graph),
                    Some(projection),
                    &generation,
                    &UVForestExportSettings { computed: true },
                    |_, atom| atom,
                )?;
                let node_indices = exported
                    .node_terms
                    .iter()
                    .map(|term| term.node_index)
                    .collect::<std::collections::BTreeSet<_>>();
                // Require the empty-forest root and at least one actual UV node;
                // merely exporting the uncomputed forest DOT misses this failure.
                assert!(
                    node_indices.len() > 1,
                    "{} (summed={explicit_sum}, projected={projected}): no nonempty computed forest",
                    generation.uv.orchestrator
                );
                assert!(exported.node_terms.iter().all(|term| {
                    term.forest_index == 0
                        && term.residue_index == CutCFFIndex::new_all_none()
                        && term.dot.contains("forest_residue_index")
                }));

                let topology = export_graph(
                    &graph,
                    CutStructure::empty(&graph),
                    None,
                    &generation,
                    &UVForestExportSettings { computed: false },
                    |_, atom| atom,
                )?;
                assert!(!topology.forest_dot.is_empty());
                assert!(topology.node_terms.is_empty());
            }
        }
        Ok(())
    }
}
