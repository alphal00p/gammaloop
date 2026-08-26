use color_eyre::Result;
use eyre::{Context, eyre};
use linnet::half_edge::subgraph::SubSetOps;
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use symbolica::atom::{Atom, AtomCore};

use crate::{
    cff::CutCFFIndex,
    graph::{
        Graph,
        cuts::{CutSet, ResidueSelector},
        parse::string_utils::dot_statement_value,
    },
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
    pub(crate) fn export_uv_forest_graph(
        &self,
        graph_id: usize,
        orientation: Option<OrientationProjection<'_>>,
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
                export_graph(
                    &term.graph,
                    cut_structure,
                    orientation,
                    generation_settings,
                    export_settings,
                    |_, atom| atom,
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
                    export_settings,
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
    let cut_woods = CutWoods::new(cut_structure, graph, generation_settings);
    let mut cut_forests = cut_woods.unfold(graph);
    let Some(forest) = cut_forests.forests.first_mut() else {
        return Err(eyre!("Legacy UV exporter produced no forest"));
    };
    forest_dot.push_str(&forest.dot_serialize_for_export(forest_name));
    forest_dot.push('\n');

    if !export_settings.computed {
        return Ok(());
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
    let wood = HedgePosetWood::new(cut_structure, graph, generation_settings);
    let mut forests = wood.unfold();
    forest_dot.push_str(&name_dot_graph(forests.dot_serialize(), forest_name));
    forest_dot.push('\n');

    if !export_settings.computed {
        return Ok(());
    }

    forests.compute(
        graph,
        crate::utils::vakint()?,
        orientation.ok_or_else(|| {
            eyre!("Computed UV forest export requires its stored production CFF expression")
        })?,
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
    let full_num = term
        .numerator
        .printer(SpensoPrintSettings::typst_options())
        .to_string();
    dot_graph
        .global_data
        .statements
        .insert("full_num".into(), dot_statement_value(&full_num));

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
    use symbolica::{atom::Atom, function};

    use super::{
        UVForestExportSettings, UVForestNodeExpression, UVForestNodeTerm, export_graph,
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
            UVOrchestrator,
            approx::{CutStructure, OrientationProjection},
        },
    };

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
        test_initialise().unwrap();
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
        })
        .unwrap();
        let term = UVForestNodeExpression {
            forest_index: 0,
            node_index: 1,
            node_key: "{1}".to_string(),
            term_index: 2,
            residue_index: CutCFFIndex::new_all_none(),
            numerator: function!(GS.uv_truncate, Atom::num(1)),
        };

        let exported = node_expression_to_dot(&graph, "uv_typst", term).unwrap();

        assert!(exported.dot.contains(r#"full_num = "op(\"Tr\")(1)";"#));
        assert!(!exported.dot.contains("gammalooprs::Truncate"));
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
                let (graph, production, expected) = match &process.collection {
                    ProcessCollection::Amplitudes(amplitudes) => {
                        let graph = &amplitudes["default"].graphs[0];
                        (
                            &graph.graph,
                            graph.derived_data.cff_expression.as_ref().unwrap(),
                            graph.derived_data.all_mighty_integrand.clone(),
                        )
                    }
                    ProcessCollection::CrossSections(cross_sections) => {
                        let graph = &cross_sections["default"].supergraphs[0];
                        assert_eq!(graph.cuts.len(), 1);
                        let integrands = &graph.derived_data.cut_paramatric_integrand
                            [crate::processes::cross_section::CutGroupId(0)]
                        .integrands;
                        assert_eq!(integrands.iter().count(), 1);
                        (
                            &graph.graph,
                            graph.derived_data.global_cff_expression.as_ref().unwrap(),
                            integrands.iter().next().unwrap().1.clone(),
                        )
                    }
                };
                let options = graph.production_cff_3d_expression_options(&settings.generation)?;
                let orientation = OrientationProjection::exact_expression(
                    production,
                    &options,
                    &settings.generation.orientation_pattern,
                    true,
                );
                let exported = process
                    .get_integrand("default")?
                    .require_generated()?
                    .export_uv_forest_graph(
                        0,
                        Some(orientation),
                        &settings.generation,
                        &UVForestExportSettings { computed: true },
                    )?;
                let [root] = exported.node_terms.as_slice() else {
                    panic!("unsubtracted scalar fixture must export exactly its root");
                };
                assert_eq!(root.node_index, 0);
                let parsed = Graph::from_string(&root.dot, &model)?;
                let split = graph
                    .iter_edges_of(
                        &graph
                            .full_filter()
                            .subtract(&graph.initial_state_cut)
                            .subtract(&graph.tree_edges),
                    )
                    .map(|(_, edge, _)| GS.split_mom_pattern_simple(edge))
                    .collect::<Vec<_>>();
                let actual = UvMarker::new(&settings.generation.uv).finish(
                    &parsed[0]
                        .global_prefactor
                        .num
                        .replace_multiple(&split)
                        .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                        .with(W_.d_),
                );
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
