use std::collections::BTreeSet;

use color_eyre::eyre::{Result, WrapErr, ensure, eyre};
use three_dimensional_reps::{ThermalDistributionFactor, utils::rank_i64};

use crate::{
    graph::{Graph, LoopMomentumBasis, lmb::LMBwithEdges},
    momentum::{SignOrZero, sample::LoopIndex},
    utils::hyperdual_utils::mixed_derivative_shape,
    uv::uv_graph::UVE,
};

use super::FermiSurfaceProduct;

impl FermiSurfaceProduct {
    /// Choose a graph basis whose first loop momenta are the selected fermion
    /// momenta. Independent internal routes extend to a unimodular graph basis;
    /// its complete edge signatures retain the original external shifts.
    pub fn new(graph: &Graph, factors: &[ThermalDistributionFactor]) -> Result<Self> {
        ensure!(
            !factors.is_empty(),
            "Fermi-surface product must not be empty"
        );
        let parent = &graph.loop_momentum_basis;
        let mut edges = BTreeSet::new();
        let mut rows = Vec::with_capacity(factors.len());
        for factor in factors {
            ensure!(
                factor.derivative_order > 0,
                "Fermi-surface edge {} requires a positive distribution derivative order",
                factor.edge_id
            );
            ensure!(
                matches!(factor.sign, -1 | 1),
                "Fermi-surface edge {} has invalid thermal sign {}",
                factor.edge_id,
                factor.sign
            );
            ensure!(
                edges.insert(factor.edge_id),
                "Fermi-surface product repeats edge {}",
                factor.edge_id
            );
            let (pair, _, edge) = graph
                .iter_edges()
                .find(|(_, edge_id, _)| *edge_id == factor.edge_id)
                .ok_or_else(|| eyre!("Unknown Fermi-surface edge {}", factor.edge_id))?;
            ensure!(
                pair.is_paired() && edge.data.is_fermion(),
                "Fermi-surface edge {} must be an internal fermion edge",
                factor.edge_id
            );
            ensure!(
                edge.data.chemical_potential_atom().is_some(),
                "Fermi-surface edge {} has no chemical potential",
                factor.edge_id
            );
            let signature = parent.edge_signatures.get(factor.edge_id).ok_or_else(|| {
                eyre!("Missing routing for Fermi-surface edge {}", factor.edge_id)
            })?;
            ensure!(
                signature.internal.len() == parent.loop_edges.len()
                    && signature.external.len() == parent.ext_edges.len(),
                "Incomplete routing for Fermi-surface edge {}",
                factor.edge_id
            );
            rows.push(
                signature
                    .internal
                    .iter()
                    .map(|sign| *sign * 1_i64)
                    .collect(),
            );
        }
        // Fixed external shifts cannot make dependent integration directions
        // independent. Check before asking the spanning-forest builder to omit
        // these edges, because a dependent choice can disconnect its guide.
        ensure!(
            rank_i64(&rows) == factors.len(),
            "Fermi-surface product requires independent loop-momentum routes; selected edges {:?} are dependent",
            edges
        );
        let selected_edges = factors
            .iter()
            .map(|factor| factor.edge_id)
            .collect::<Vec<_>>();
        // Complete the active directions with independent fermion routes
        // before choosing arbitrary spanning-tree chords. Otherwise an
        // untouched occupation step can acquire an artificial dependence on
        // a localized radius (e.g. the other fermion of a dotted sunset).
        let mut basis_edges = selected_edges.clone();
        for (pair, edge_id, edge) in graph.iter_edges() {
            if basis_edges.len() == parent.loop_edges.len() {
                break;
            }
            if !pair.is_paired() || !edge.data.is_fermion() {
                continue;
            }
            rows.push(
                parent.edge_signatures[edge_id]
                    .internal
                    .iter()
                    .map(|sign| *sign * 1_i64)
                    .collect(),
            );
            if rank_i64(&rows) == rows.len() {
                basis_edges.push(edge_id);
            } else {
                rows.pop();
            }
        }
        let mut lmb = graph
            .lmb_with_loop_edges(basis_edges.as_slice())
            .wrap_err("Could not construct an independent Fermi-surface basis")?;
        graph.canonicalize_lmb_external_order(&mut lmb);
        ensure!(
            lmb.loop_edges.len() == parent.loop_edges.len() && lmb.ext_edges == parent.ext_edges,
            "Fermi-surface basis does not preserve the graph's loop and external momentum spaces"
        );
        for (slot, edge) in selected_edges.iter().enumerate() {
            let position = lmb
                .loop_edges
                .iter()
                .position(|candidate| candidate == edge)
                .ok_or_else(|| {
                    eyre!("Fermi-surface edge {edge} is absent from the adapted basis")
                })?;
            lmb.swap_loops(LoopIndex(slot), LoopIndex(position));
        }
        for (slot, edge) in selected_edges.iter().enumerate() {
            let signature = &lmb.edge_signatures[*edge];
            ensure!(
                signature.internal.iter().enumerate().all(|(index, sign)| {
                    *sign
                        == if index == slot {
                            SignOrZero::Plus
                        } else {
                            SignOrZero::Zero
                        }
                }) && signature
                    .external
                    .iter()
                    .all(|sign| *sign == SignOrZero::Zero),
                "Fermi-surface edge {edge} is not a defining momentum in its adapted basis"
            );
        }
        let derivative_orders = factors
            .iter()
            .map(|factor| factor.derivative_order - 1)
            .collect::<Vec<_>>();
        let shape = mixed_derivative_shape(&derivative_orders);
        Ok(Self {
            factors: factors.to_vec(),
            lmb,
            derivative_orders,
            shape,
        })
    }

    pub fn lmb(&self) -> &LoopMomentumBasis {
        &self.lmb
    }

    pub fn factors(&self) -> &[ThermalDistributionFactor] {
        &self.factors
    }
}

#[cfg(test)]
pub(super) fn test_graph() -> Result<Graph> {
    use crate::{dot, graph::parse::IntoGraph, initialisation::test_initialise};

    test_initialise()?;
    dot!(digraph shifted_fermi_cycles {
        node [num=1]
        edge [num=1 particle="d"]
        ext [style=invis]
        A -> B [id=0 lmb_id=0]
        B -> A [id=1]
        A -> C [id=2 lmb_id=1]
        C -> A [id=3]
        ext -> B [id=4 particle="H"]
        C -> ext [id=5 particle="H"]
    })
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::involution::EdgeIndex;

    use super::*;
    use crate::{
        dot, graph::parse::IntoGraph, initialisation::test_initialise, momentum::ThreeMomentum,
        utils::F,
    };

    #[test]
    fn independent_fermi_basis_preserves_affine_routes() -> Result<()> {
        let graph = test_graph()?;
        let factors = [
            ThermalDistributionFactor {
                edge_id: EdgeIndex(3),
                sign: -1,
                derivative_order: 2,
            },
            ThermalDistributionFactor {
                edge_id: EdgeIndex(1),
                sign: 1,
                derivative_order: 3,
            },
        ];
        let product = FermiSurfaceProduct::new(&graph, &factors)?;
        assert_eq!(product.factors(), factors);
        assert_eq!(product.lmb().loop_edges.raw, [EdgeIndex(3), EdgeIndex(1)]);
        assert_eq!(product.derivative_orders, [1, 2]);
        assert_eq!(product.shape.len(), 6);
        let parent = &graph.loop_momentum_basis;
        assert!(factors.iter().any(|factor| {
            parent.edge_signatures[factor.edge_id]
                .external
                .iter()
                .any(|sign| *sign != SignOrZero::Zero)
        }));

        let loops = [
            ThreeMomentum::new(F(1.0), F(2.0), F(3.0)),
            ThreeMomentum::new(F(4.0), F(5.0), F(6.0)),
        ];
        let externals = vec![ThreeMomentum::new(F(2.0), F(-3.0), F(4.0)); parent.ext_edges.len()];
        let adapted = product
            .lmb()
            .loop_edges
            .iter()
            .map(|edge| parent.edge_signatures[*edge].compute_momentum_untyped(&loops, &externals))
            .collect::<Vec<_>>();
        for (edge, signature) in &parent.edge_signatures {
            assert_eq!(
                signature.compute_momentum_untyped(&loops, &externals),
                product.lmb().edge_signatures[edge].compute_momentum_untyped(&adapted, &externals),
                "affine route for edge {edge}"
            );
        }

        let partial = FermiSurfaceProduct::new(&graph, &factors[..1])?;
        assert_eq!(partial.lmb().loop_edges.len(), 2);
        assert_eq!(partial.lmb().loop_edges[LoopIndex(0)], EdgeIndex(3));
        Ok(())
    }

    #[test]
    fn fermi_basis_keeps_the_other_sunset_occupation_independent() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph nonaligned_dotted_sunset {
            node [num=1]
            edge [num=1]
            A -> C [id=0 particle="b"]
            C -> B [id=1 particle="b"]
            B -> A [id=2 particle="b" lmb_id=0]
            A -> D [id=3 particle="H" lmb_id=1]
            D -> B [id=4 particle="H"]
        })?;
        let parent = &graph.loop_momentum_basis;
        for edge in [EdgeIndex(0), EdgeIndex(1)] {
            assert!(!parent.loop_edges.contains(&edge));
            assert!(
                parent.edge_signatures[edge]
                    .internal
                    .iter()
                    .all(|sign| *sign != SignOrZero::Zero)
            );
            let product = FermiSurfaceProduct::new(
                &graph,
                &[ThermalDistributionFactor {
                    edge_id: edge,
                    sign: 1,
                    derivative_order: 1,
                }],
            )?;
            assert_eq!(product.lmb().loop_edges[LoopIndex(0)], edge);
            assert_eq!(
                product.lmb().edge_signatures[EdgeIndex(2)].internal[LoopIndex(0)],
                SignOrZero::Zero,
                "the other occupation step must not depend on the localized radius"
            );
        }
        Ok(())
    }

    #[test]
    fn fermi_product_rejects_dependent_routes_despite_external_shifts() -> Result<()> {
        let graph = test_graph()?;
        let parent = &graph.loop_momentum_basis;
        assert_eq!(
            parent.edge_signatures[EdgeIndex(0)].internal,
            parent.edge_signatures[EdgeIndex(1)].internal
        );
        assert_ne!(
            parent.edge_signatures[EdgeIndex(0)].external,
            parent.edge_signatures[EdgeIndex(1)].external
        );
        let factors = [0, 1].map(|edge| ThermalDistributionFactor {
            edge_id: EdgeIndex(edge),
            sign: 1,
            derivative_order: 1,
        });
        let error = FermiSurfaceProduct::new(&graph, &factors).err().unwrap();
        assert!(error.to_string().contains("dependent"));
        let duplicate = [factors[0], factors[0]];
        let error = FermiSurfaceProduct::new(&graph, &duplicate).err().unwrap();
        assert!(error.to_string().contains("repeats edge"));
        Ok(())
    }

    #[test]
    fn fermi_product_rejects_malformed_factors_and_accepts_a_tadpole() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph fermi_tadpoles {
            node [num=1]
            edge [num=1]
            A -> A [id=0 particle="d"]
            A -> A [id=1 particle="g"]
        })?;
        let factor = ThermalDistributionFactor {
            edge_id: EdgeIndex(0),
            sign: 1,
            derivative_order: 1,
        };
        assert!(FermiSurfaceProduct::new(&graph, &[factor]).is_ok());
        assert!(FermiSurfaceProduct::new(&graph, &[]).is_err());
        for invalid in [
            ThermalDistributionFactor {
                derivative_order: 0,
                ..factor
            },
            ThermalDistributionFactor { sign: 0, ..factor },
            ThermalDistributionFactor {
                edge_id: EdgeIndex(1),
                ..factor
            },
            ThermalDistributionFactor {
                edge_id: EdgeIndex(20),
                ..factor
            },
        ] {
            assert!(FermiSurfaceProduct::new(&graph, &[invalid]).is_err());
        }
        Ok(())
    }
}
