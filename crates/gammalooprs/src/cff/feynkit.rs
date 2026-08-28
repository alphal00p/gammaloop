use std::collections::BTreeSet;

use color_eyre::Result;
use eyre::eyre;
use feynkit_cff::{CffOptions, CffResult, EdgeId, HedgeEdgeRole, HedgeGraphCffExt, ShiftRewrite};
use linnet::half_edge::involution::EdgeIndex;

use crate::{
    cff::generation::GammaLoopGraphCffExt, graph::Graph, settings::global::OrientationPattern,
};

impl GammaLoopGraphCffExt for Graph {
    fn generate_cff(
        &mut self,
        contract_edges: &[EdgeIndex],
        shift_rewrite: &Option<ShiftRewrite>,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CffResult> {
        let initial_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge, _)| edge)
            .collect::<BTreeSet<_>>();
        let contracted_edges = contract_edges.iter().copied().collect::<BTreeSet<_>>();
        if let Some(edge) = contracted_edges.intersection(&initial_edges).next() {
            return Err(eyre!(
                "initial-state edge {edge} cannot be contracted during root CFF generation"
            ));
        }

        let dummy_edges = self
            .iter_edges()
            .filter_map(|(_, edge, data)| data.data.is_dummy.then_some(edge))
            .collect::<BTreeSet<_>>();
        let external_edges = self
            .iter_edges()
            .filter_map(|(pair, edge, _)| pair.is_unpaired().then_some(edge))
            .collect::<BTreeSet<_>>();

        let mut options = contract_edges
            .iter()
            .fold(CffOptions::default(), |options, edge| {
                options.with_contracted_edge(EdgeId::new(edge.0))
            });
        if let Some(rewrite) = shift_rewrite {
            options = options.with_shift_rewrite(rewrite.clone());
        }

        let full_graph = self.underlying.full_filter();
        let mut result = self.underlying.build_cff_from_subgraph_with_edge_roles(
            &full_graph,
            options,
            &[],
            |edge| {
                if dummy_edges.contains(&edge) {
                    HedgeEdgeRole::Omitted
                } else if initial_edges.contains(&edge) {
                    HedgeEdgeRole::InitialState
                } else if external_edges.contains(&edge) {
                    HedgeEdgeRole::UnorientedExternal
                } else {
                    HedgeEdgeRole::Standard
                }
            },
            &mut self.surface_cache,
        )?;
        result
            .expression
            .retain_orientations(|orientation| orientation_pattern.filter(orientation));
        result.report.acyclic_orientations = result.expression.orientations().len();
        result.report.unfolded_terms = result.expression.unfolded_term_count();
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use feynkit_cff::{
        CffError, CffOptions, EnergySurface, HSurface, HedgeGraphCffExt, ShiftRewrite, Surface,
        VertexSet,
    };
    use feynkit_generator::{GenerationOptions, Generator, Process};
    use linnet::half_edge::{
        involution::EdgeIndex,
        subgraph::{ModifySubSet, SuBitGraph},
    };

    use crate::{
        cff::generation::GammaLoopGraphCffExt,
        feyngen::feynkit::GenerationResultGammaLoopExt,
        finalized_runtime_dot,
        graph::{Graph, parse::IntoFinalizedRuntimeGraph},
        initialisation::test_initialise,
        settings::global::OrientationPattern,
        utils::load_generic_model,
    };

    fn bubble() -> Graph {
        finalized_runtime_dot!(
            digraph bubble {
                projector=1
                edge [particle=scalar_1 num=1]
                node [num=1]
                ext [style=invis]
                ext -> A:0 [id=3 sink="{ufo_order:0}"]
                B:1 -> ext [id=2 source="{ufo_order:2}"]
                A:2 -> B:3 [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:0}"]
                A:4 -> B:5 [id=0 source="{ufo_order:2}" sink="{ufo_order:1}"]
            },
            "scalars"
        )
        .unwrap()
    }

    fn cross_section_bubble() -> Graph {
        let model = Arc::new(load_generic_model("scalars"));
        let process = Process::cross_section(["scalar_1"], ["scalar_1", "scalar_1"])
            .with_loop_count(1, 1)
            .unwrap();
        let generated = Generator::new(model)
            .generate(
                &process,
                &GenerationOptions::default().threads(1).max_vertices(2),
            )
            .unwrap();
        generated
            .to_gamma_loop_graphs()
            .unwrap()
            .into_iter()
            .next()
            .expect("the scalar one-loop cross section must generate a runtime graph")
    }

    fn triangle() -> Graph {
        finalized_runtime_dot!(
            digraph triangle {
                projector=1
                edge [particle=scalar_1 num=1]
                node [num=1]
                ext [style=invis]
                ext -> A:0 [id=4 sink="{ufo_order:0}"]
                C:1 -> ext [id=3 source="{ufo_order:2}"]
                A:2 -> B:3 [id=2 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:0}"]
                B:4 -> C:5 [id=1 source="{ufo_order:1}" sink="{ufo_order:0}"]
                A:6 -> C:7 [id=0 source="{ufo_order:2}" sink="{ufo_order:1}"]
            },
            "scalars"
        )
        .unwrap()
    }

    fn triangle_left_subgraph(graph: &Graph) -> SuBitGraph {
        let mut subgraph: SuBitGraph = graph.empty_subgraph();
        for (_, neighbours, vertex) in graph.iter_nodes() {
            if vertex.name.as_str() == "A" || vertex.name.as_str() == "B" {
                neighbours.for_each(|hedge| subgraph.add(hedge));
            }
        }
        subgraph
    }

    #[test]
    fn root_amplitude_default_snapshot() {
        test_initialise().unwrap();
        let mut graph = bubble();
        let expression = graph
            .generate_cff(&[], &None, &OrientationPattern::default())
            .unwrap();

        insta::assert_json_snapshot!(
            "root_amplitude_default",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_amplitude_selected_orientation_snapshot() {
        test_initialise().unwrap();
        // FeynKit enumerates internal-edge orientations from a bitset: unset bits retain
        // the graph direction and set bits reverse it. Fixed external and boundary directions
        // pad the GammaLoop orientation, and cyclic orientations are filtered before conversion.
        let selected = OrientationPattern::from_user_pattern("(1,1,0,0)").unwrap();
        let mut graph = bubble();
        let expression = graph.generate_cff(&[], &None, &selected).unwrap();
        assert_eq!(expression.orientations.len(), 1);

        insta::assert_json_snapshot!(
            "root_amplitude_selected_orientation",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_tree_snapshot() {
        test_initialise().unwrap();
        let mut graph = triangle();
        let expression = graph
            .generate_cff(&[], &None, &OrientationPattern::default())
            .unwrap();

        insta::assert_json_snapshot!(
            "root_tree",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_contracted_edge_snapshot() {
        test_initialise().unwrap();
        let mut graph = triangle();
        let expression = graph
            .generate_cff(&[EdgeIndex(1)], &None, &OrientationPattern::default())
            .unwrap();

        insta::assert_json_snapshot!(
            "root_contracted_edge",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_initial_state_shift_and_cache_snapshot() {
        test_initialise().unwrap();
        let mut graph = cross_section_bubble();
        graph.surface_cache.intern(Surface::Energy(EnergySurface {
            energies: vec![EdgeIndex(99)],
            external_shift: vec![].into(),
            vertex_set: VertexSet::dummy(),
        }));
        graph.surface_cache.intern(Surface::H(HSurface {
            positive_energies: vec![EdgeIndex(98)],
            negative_energies: vec![EdgeIndex(99)],
            external_shift: vec![].into(),
            vertex_set: VertexSet::dummy(),
        }));
        let initial_edge = graph
            .iter_edges_of(&graph.initial_state_cut)
            .next()
            .unwrap()
            .1;
        let shift_rewrite = Some(ShiftRewrite {
            dependent_momentum: initial_edge,
            replacement: vec![].into(),
        });

        let expression = graph
            .generate_cff(&[], &shift_rewrite, &OrientationPattern::default())
            .unwrap();
        insta::assert_json_snapshot!(
            "root_initial_state_shift_and_cache",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_generation_rejects_contracted_initial_state_edges() {
        test_initialise().unwrap();
        let mut graph = cross_section_bubble();
        let initial_edge = graph
            .iter_edges_of(&graph.initial_state_cut)
            .next()
            .unwrap()
            .1;

        let error = graph
            .generate_cff(&[initial_edge], &None, &OrientationPattern::default())
            .unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "initial-state edge {initial_edge} cannot be contracted during root CFF generation"
            )
        );
    }

    #[test]
    fn subgraph_reversed_boundary_snapshot() {
        test_initialise().unwrap();
        let graph = triangle();
        let subgraph = triangle_left_subgraph(&graph);
        let reversed_dangling = [EdgeIndex(0)];
        let initial_edges = graph.get_edges_in_initial_state_cut();
        let options = initial_edges
            .into_iter()
            .fold(CffOptions::default(), |options, edge| {
                options.with_initial_state_edge(edge)
            });
        let mut surface_cache = graph.surface_cache.clone();
        let expression = graph
            .underlying
            .build_cff_from_subgraph(&subgraph, options, &reversed_dangling, &mut surface_cache)
            .unwrap();

        insta::assert_json_snapshot!(
            "subgraph_reversed_boundary",
            serde_json::json!({
                "expression": expression,
                "surface_cache": surface_cache,
            })
        );
    }

    #[test]
    fn empty_subgraph_reports_the_typed_feynkit_error() {
        test_initialise().unwrap();
        let graph = triangle();
        let empty: SuBitGraph = graph.empty_subgraph();
        let mut cache = graph.surface_cache.clone();

        let error = graph
            .underlying
            .build_cff_from_subgraph(&empty, CffOptions::default(), &[], &mut cache)
            .unwrap_err();

        assert_eq!(error, CffError::EmptyGraph);
    }
}
