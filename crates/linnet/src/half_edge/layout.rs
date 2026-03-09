use super::{nodestore::NodeStorageVec, HedgeGraph};

pub mod spring;

pub mod simulatedanneale;

pub mod force;

impl<E, V, H> HedgeGraph<E, V, H, NodeStorageVec<V>> {
    // pub fn layout(
    //     self,
    //     mut settings: LayoutSettings,
    // ) -> HedgeGraph<LayoutEdge<E>, LayoutVertex<V>, H, NodeStorageVec<LayoutVertex<V>>> {
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
pub mod test {}
