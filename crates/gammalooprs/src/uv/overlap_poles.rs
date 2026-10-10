//! Physical, signed FullH input for an explicit offline soft-overlap proof.
//! Never resolve retained graph numerators into the complete forest expression.

use super::*;
use crate::{
    graph::FeynmanGraph,
    uv::{Integrands, approx::direct_3d::physical_full_h_fixture},
};

#[derive(Clone, Debug)]
pub(super) struct OverlapStep {
    /// Stored Foata level, in dependency order from the innermost operation.
    pub(super) level: usize,
    pub(super) current_edges: Vec<usize>,
    pub(super) given_edges: Vec<usize>,
    pub(super) scheme: ApproximationType,
    pub(super) dod: i32,
}

pub(super) struct FinalizedForestNode {
    /// Actual unfolded graph identity, not its position in `nodes`.
    pub(super) node_index: usize,
    pub(super) node_key: String,
    pub(super) steps: Vec<OverlapStep>,
    /// Documentary only: all subtraction signs are already in `integrands`.
    pub(super) forest_sign: i32,
    /// Complete orientation sum with literal, separately retained FnMapEntry bodies.
    pub(super) integrands: Integrands,
}

pub(super) struct PhysicalOverlapInput {
    pub(super) graph: Graph,
    pub(super) nodes: Vec<FinalizedForestNode>,
    pub(super) orientation_count: usize,
    pub(super) negative_control_node: usize,
}

impl PhysicalOverlapInput {
    /// Generate symbolic FullH only. No tensor evaluator is constructed here.
    pub(super) fn generate() -> Result<Self> {
        let (mut graph, generation) = physical_full_h_fixture()?;
        assert!(generation.explicit_orientation_sum_only);
        assert!(generation.orientation_pattern.pat.is_none());
        let options = graph.production_cff_3d_expression_options(&generation)?;
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let contract_edges = graph
            .iter_edges_of(
                &graph
                    .tree_edges
                    .subtract(&graph.initial_state_cut)
                    .subtract(&graph.external_filter::<SuBitGraph>()),
            )
            .map(|(_, edge, _)| edge)
            .collect_vec();
        let production = graph.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonization,
            &options,
            None,
        )?;
        let orientation_count = production.expression.orientations.len();
        assert_eq!(orientation_count, 102);
        debug_tags!(#uv, #overlap_poles;
            stage = "overlap_poles_orientation_inventory", graph = %graph.name,
            orientation_count,
            file.orientations = ?production.expression.orientations.iter().map(|row| {
                (&row.data, &row.loop_energy_map, &row.edge_energy_map)
            }).collect_vec(),
            "Complete unfiltered production orientation and energy-map inventory"
        );
        let orientation = OrientationProjection::exact_expression(
            &production,
            &options,
            &generation.orientation_pattern,
            true,
        );
        let mut forests = Wood::new(CutStructure::empty(&graph), &graph, &generation.uv)?.unfold();
        let [(compatible, cutset)] = forests.cuts.as_slice() else {
            return Err(eyre!("the overlap fixture must have one empty cut"));
        };
        let cutset = cutset.clone();
        let order = forests.compatible_topological_order(compatible)?;
        let mut histogram = BTreeMap::<usize, usize>::new();
        let mut inventory = Vec::new();
        for &node in &order {
            let operation = &forests.graph[node];
            let steps = operation
                .key
                .iter_levels_top_down()
                .enumerate()
                .flat_map(|(level, operations)| {
                    operations.iter_leaf_ops().map(move |op| (level, op))
                })
                .map(|(level, op)| {
                    let (current, given) = forests.wood.current_given_pair(op.data, level);
                    let edges = |subgraph: &SuBitGraph| {
                        graph
                            .iter_edges_of(subgraph)
                            .map(|(_, edge, _)| usize::from(edge))
                            .sorted()
                            .collect_vec()
                    };
                    OverlapStep {
                        level,
                        current_edges: edges(current.subgraph()),
                        given_edges: edges(given.subgraph()),
                        scheme: current.renormalization_scheme(),
                        dod: current.dod(),
                    }
                })
                .collect_vec();
            assert_eq!(steps.len(), operation.key.op_count());
            *histogram.entry(steps.len()).or_default() += 1;
            let forest_sign = if steps.len().is_multiple_of(2) { 1 } else { -1 };
            debug_tags!(#uv, #overlap_poles;
                stage = "overlap_poles_node_inventory", graph = %graph.name,
                node_index = usize::from(node), node_key = %operation,
                foata_levels = %operation.foata_level_labels(),
                steps = ?steps.iter().map(|step| (
                    step.level, &step.current_edges, &step.given_edges, step.scheme, step.dod,
                )).collect_vec(),
                forest_sign, sign_already_in_expression = true,
                "Compatible signed forest node before symbolic generation"
            );
            inventory.push((node, steps, forest_sign));
        }
        assert_eq!(histogram.get(&0), Some(&1));
        // Select by actual nested subgraphs/schemes, never a transient node number.
        let control_chain = [vec![8, 9], vec![2, 5, 6, 7, 8, 9], (2..10).collect_vec()];
        let controls = inventory
            .iter()
            .filter(|(_, steps, _)| {
                steps.len() == 3
                    && steps
                        .iter()
                        .zip(&control_chain)
                        .all(|(step, edges)| &step.current_edges == edges)
                    && steps.iter().map(|step| step.scheme).collect_vec()
                        == [
                            ApproximationType::IR,
                            ApproximationType::MUV,
                            ApproximationType::IR,
                        ]
                    && steps.iter().map(|step| step.dod).collect_vec() == [2, 0, 2]
            })
            .collect_vec();
        let [(negative_control, _, sign)] = controls.as_slice() else {
            return Err(eyre!(
                "expected exactly one H2/U0/H2 triple-chain control, found {}",
                controls.len()
            ));
        };
        assert_eq!(*sign, -1);
        let negative_control_node = usize::from(*negative_control);
        debug_tags!(#uv, #overlap_poles;
            stage = "overlap_poles_generation_start", graph = %graph.name,
            orientation_count, forest_nodes = inventory.len(), operation_count_histogram = ?histogram,
            negative_control_node, negative_control_chain = ?control_chain,
            file.generation_settings = ?generation,
            "Generate the complete FullH signed forest; retain all numerator definitions"
        );
        forests.compute(
            &mut graph,
            crate::utils::vakint()?,
            orientation,
            &generation.uv,
        )?;
        let mut nodes = Vec::with_capacity(inventory.len());
        for (node, steps, forest_sign) in inventory {
            let operation = &forests.graph[node];
            let computed = forests
                .compute_store
                .require(operation)?
                .cut(operation, &cutset)?;
            // Match production's per-node color collection before the forest sum.
            // FinalIntegrands::iter/export would resolve every numerator; do not use them.
            let integrands = computed
                .final_integrands
                .map_expressions(|atom| Ok(atom.collect_color()))?
                .into_integrands();
            assert_eq!(integrands.iter().count(), 1);
            assert!(
                integrands
                    .iter()
                    .all(|(key, _)| *key == crate::cff::CutCFFIndex::new_all_none())
            );
            let provenance = forests.node_export_3d_provenance(node, &computed.local_3d);
            debug_tags!(#uv, #overlap_poles;
                stage = "overlap_poles_node_captured", node_index = usize::from(node),
                node_key = %operation, operation_count = steps.len(), forest_sign,
                root_bytes = integrands.iter().map(|(_, atom)| atom.as_view().get_byte_size()).sum::<usize>(),
                retained_definitions = integrands.numerators().len(),
                file.provenance = ?provenance,
                "Captured finalized roots and separate literal numerator definitions"
            );
            nodes.push(FinalizedForestNode {
                node_index: usize::from(node),
                node_key: operation.to_string(),
                steps,
                forest_sign,
                integrands,
            });
        }
        assert!(
            nodes
                .iter()
                .any(|node| !node.integrands.numerators().is_empty())
        );
        assert_eq!(nodes.iter().filter(|node| node.steps.is_empty()).count(), 1);
        Ok(Self {
            graph,
            nodes,
            orientation_count,
            negative_control_node,
        })
    }

    /// Match the production post-sum momentum split for the full sum or the
    /// explicit omit-one-node control. Retained RHS remain separate throughout.
    pub(super) fn sum_except(&self, omit_node: Option<usize>) -> Result<Integrands> {
        if let Some(omit) = omit_node {
            eyre::ensure!(
                self.nodes.iter().any(|node| node.node_index == omit),
                "unknown omitted forest node {omit}"
            );
        }
        let mut nodes = self
            .nodes
            .iter()
            .filter(|node| Some(node.node_index) != omit_node);
        let first = nodes
            .next()
            .ok_or_else(|| eyre!("no forest nodes remain"))?;
        let sum = first
            .integrands
            .clone()
            .zip_add(nodes.map(|node| node.integrands.clone()))?;
        self.finish(sum)
    }

    /// Diagnostic contribution in the same terminal momentum convention as the sum.
    /// Its subtraction sign is already included; definitions remain literal.
    pub(super) fn node_integrands(&self, node_index: usize) -> Result<Integrands> {
        let node = self
            .nodes
            .iter()
            .find(|node| node.node_index == node_index)
            .ok_or_else(|| eyre!("unknown forest node {node_index}"))?;
        self.finish(node.integrands.clone())
    }

    fn finish(&self, integrands: Integrands) -> Result<Integrands> {
        let split = self
            .graph
            .iter_edges_of(
                &self
                    .graph
                    .full_filter()
                    .subtract(&self.graph.initial_state_cut)
                    .subtract(&self.graph.tree_edges),
            )
            .filter_map(|(pair, edge, _)| {
                pair.is_paired().then(|| GS.split_mom_pattern_simple(edge))
            })
            .collect_vec();
        let finish = |atom: &Atom| {
            Ok(atom
                .replace_multiple(&split)
                .replace(function!(GS.den, W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_))
        };
        integrands.fallible_map(finish)?.map_numerators(finish)
    }
}

#[path = "overlap_proof.rs"]
mod proof;

#[test]
#[ignore = "offline physical FullH soft-overlap proof; generates all 102 orientations without evaluators"]
fn physical_full_h_overlap_poles() -> Result<()> {
    let input = PhysicalOverlapInput::generate()?;
    proof::check(&input)
}
