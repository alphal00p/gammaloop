use clarabel::algebra::*;
use clarabel::solver::*;
use lorentz_vector::LorentzVector;

use crate::graph::Graph;

pub fn construct_problem(
    graph: &Graph,
    esurface_ids: &[usize],
    external_momenta: &[LorentzVector<f64>],
) -> DefaultSolver {
    let lmb = graph.loop_momentum_basis;
    let num_loops = lmb.basis.len();
    let num_loop_vars = 3 * num_loops;

    let esurfaces = &graph
        .derived_data
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    let esurface_derived_data = graph.derived_data.esurface_derived_data.as_ref().unwrap();

    // upper bound on the number of propagators
    let mut inequivalent_propagators: Vec<usize> = Vec::with_capacity(lmb.edge_signatures.len());

    // for now, each propagator has it's own mass
    let mut masses = Vec::with_capacity(lmb.edge_signatures.len());

    // loop_signatures
    let mut propagator_loop_signatures: Vec<Vec<isize>> =
        Vec::with_capacity(lmb.edge_signatures.len());

    let mut num_massive_propagators = 0;

    for esurface_id in esurface_ids {
        let esurface = &esurfaces[*esurface_id];
        for &edge_id in &esurface.energies {
            if let Some(edge_position) = inequivalent_propagators
                .iter()
                .position(|&index| index == edge_id)
            {
            } else {
                // store the propagator and it's mass
                let mass = graph.edges[edge_id].particle.mass.value.map(|m| m.re);

                if mass.is_some() {
                    num_massive_propagators += 1;
                }

                inequivalent_propagators.push(edge_id);
                let loop_signature = lmb.edge_signatures[edge_id].0;
                propagator_loop_signatures.push(loop_signature);

                masses.push(mass);
            };
        }
    }

    // this is a rather ungodly expression.
    let num_constraints = 1
        + esurface_ids.len()
        + inequivalent_propagators
            .iter()
            .enumerate()
            .map(|(index, edge_id)| {
                let mut cone_size = 0;
                if masses[index].is_some() {
                    cone_size += 1;
                }

                cone_size += 3 * propagator_loop_signatures[index]
                    .iter()
                    .filter(|&x| *x != 0)
                    .count();

                cone_size
            })
            .sum::<usize>()
        + num_massive_propagators;

    todo!()
}
