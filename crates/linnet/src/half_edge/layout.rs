use super::HedgeGraph;

pub mod spring;

pub mod simulatedanneale;

pub mod force;

pub mod impred;

pub mod layered;

impl<E, V, H> HedgeGraph<E, V, H> {
    // pub fn layout(
    //     self,
    //     mut settings: LayoutSettings,
    // ) -> HedgeGraph<LayoutEdge<E>, LayoutVertex<V>, H> {
    //     let energy = GraphLayoutEnergy {
    //         graph: &self,
    //         positions: &settings.positions,
    //         params: &settings.params,
    //     };

    //     // Epoch-based annealing with adaptive delta
    //     let mut delta = settings.iters.delta; // e.g. 0.1*L
    //     let mut temp = settings.iters.temp; // e.g. 1.0
    //     let alpha_t = 0.85; // temperature decay per epoch
    //     let alpha_d = 0.7; // delta decay when acceptance too low

    //     let total = settings.iters.n_iters as usize;
    //     let epoch = 200.min(total / 5); // steps per epoch, but not too small
    //     let mut best_state = LayoutState {
    //         params: settings.init_params,
    //         delta,
    //     };
    //     let mut best_e = energy.cost(&best_state);

    //     let rng = StdRng::seed_from_u64(settings.iters.seed);
    //     let mut state = best_state.clone();

    //     for _ in (0..total).step_by(epoch) {
    //         let schedule = GeometricSchedule::new(temp, 0.95);
    //         state.delta = delta;
    //         let energy_clone = GraphLayoutEnergy {
    //             graph: energy.graph,
    //             positions: energy.positions,
    //             params: energy.params,
    //         };
    //         let mut annealer =
    //             Annealer::new(state.clone(), energy_clone, schedule, rng.clone(), epoch);
    //         let (s, e) = annealer.run();

    //         // Update best if improved
    //         if e < best_e {
    //             best_e = e;
    //             best_state = s.clone();
    //         }

    //         // Adapt step size (simplified acceptance proxy)
    //         // In practice, you'd track actual acceptance ratio from annealer
    //         let energy_improvement = (best_e - e).abs() / best_e.abs().max(1e-10);
    //         if energy_improvement < 1e-6 {
    //             delta *= alpha_d;
    //         }

    //         temp *= alpha_t;
    //         state = s;

    //         // Early termination if delta becomes too small
    //         if delta < 1e-8 {
    //             break;
    //         }
    //     }

    //     let (best_state, _best_energy) = (best_state, best_e);
    //     let best = best_state.params;

    //     // Check if we have any fixed positions
    //     let has_pins = settings
    //         .positions
    //         .vertex_positions
    //         .iter()
    //         .any(|constraints| {
    //             matches!(constraints.x, CoordinateConstraint::Fixed(_))
    //                 || matches!(constraints.y, CoordinateConstraint::Fixed(_))
    //         })
    //         || settings
    //             .positions
    //             .edge_positions
    //             .iter()
    //             .any(|(_, constraints)| {
    //                 matches!(constraints.x, CoordinateConstraint::Fixed(_))
    //                     || matches!(constraints.y, CoordinateConstraint::Fixed(_))
    //             });

    //     if has_pins {
    //         // When fixed coordinates are present, they establish the coordinate system scale
    //         // Don't rescale - keep fixed coordinates at their absolute values
    //         settings.positions.to_graph(self, &best)
    //     } else {
    //         // When no fixed coordinates are present, normalize to prevent coordinates from becoming too large
    //         let max = settings.positions.max(&best).unwrap_or(1.0);
    //         let best: Vec<_> = best.into_iter().map(|a| a / max).collect();
    //         settings.positions.scale(max);
    //         settings.positions.to_graph(self, &best)
    //     }
    // }
}

#[cfg(test)]
pub(crate) mod test {
    use cgmath::Point2;

    use super::spring::{LayoutState, PointConstraint};
    use crate::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{EdgeVec, Flow},
        nodestore::NodeForest,
        remap::NodeRemap,
        HedgeGraph, NoData, NodeIndex, NodeVec,
    };

    type Graph = HedgeGraph<PointConstraint, PointConstraint>;

    /// A five-node ring with two external legs whose nodes 1 and 3 are
    /// identified, as is and compacted. Node 3 is absorbed and the compaction
    /// moves node 4 into its slot, so the live nodes keep their order: layouts
    /// that break ties or draw random numbers in node order lay out both alike.
    pub(crate) struct PendingIdentification {
        pub(crate) graph: Graph,
        pub(crate) compacted: Graph,
        pub(crate) remap: NodeRemap,
        pub(crate) absorbed: NodeIndex,
    }

    impl PendingIdentification {
        pub(crate) fn new() -> Self {
            let mut builder = HedgeGraphBuilder::new();
            let nodes: Vec<_> = (0..5)
                .map(|_| builder.add_node(PointConstraint::default()))
                .collect();
            for (i, &node) in nodes.iter().enumerate() {
                let next = nodes[(i + 1) % nodes.len()];
                builder.add_edge(node, next, PointConstraint::default(), false);
            }
            builder.add_external_edge(nodes[0], PointConstraint::default(), false, Flow::Sink);
            builder.add_external_edge(nodes[2], PointConstraint::default(), false, Flow::Source);
            let mut graph: Graph = builder.build();
            let survivor = graph.identify_nodes(&[nodes[1], nodes[3]], PointConstraint::default());
            assert_eq!(survivor, nodes[1], "ties keep the first requested node");
            let mut compacted = graph.clone();
            let remap = compacted.compact();
            PendingIdentification {
                graph,
                compacted,
                remap,
                absorbed: nodes[3],
            }
        }

        /// The values of the live nodes, moved to their compacted slots.
        pub(crate) fn compact<T: Clone>(&self, values: &NodeVec<T>) -> NodeVec<T> {
            let mut values = values.clone();
            self.remap.apply(&mut values.0);
            values
        }

        /// Layout states of both graphs from the same distinct points; the
        /// absorbed node has a point of its own, which must not take part.
        pub(crate) fn states(
            &self,
        ) -> [LayoutState<'_, PointConstraint, PointConstraint, NoData, NodeForest>; 2] {
            let point = |t: f64| Point2::new(2.0 * t.cos() + 0.1 * t, 1.5 * t.sin());
            let nodes: NodeVec<_> = (0..self.graph.node_bound())
                .map(|i| point(i as f64))
                .collect();
            let edges: EdgeVec<_> = (0..self.graph.n_edges())
                .map(|i| point(0.7 * i as f64 + 0.3))
                .collect();
            [
                (self.graph).new_layout_state(nodes.clone(), edges.clone(), 0.5, 0.0, false),
                (self.compacted).new_layout_state(self.compact(&nodes), edges, 0.5, 0.0, false),
            ]
        }
    }
}
