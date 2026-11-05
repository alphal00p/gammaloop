//! Phase-quenched pressure differences at the completed thermal-weight boundary.

use std::collections::BTreeMap;

use linnet::half_edge::involution::{EdgeIndex, HedgePair};
use linnet::half_edge::subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike};
use symbolica::prelude::*;
use thiserror::Error;
use three_dimensional_reps::ThermalDistributionCall;
use three_dimensional_reps::thermal::{OCCUPATION_FORM, ThermalDistributionCallError};

use crate::graph::Graph;
use crate::utils::GS;
use crate::uv::uv_graph::UVE;

#[derive(Debug, Error)]
pub(crate) enum PhaseQuenchedDifferenceError {
    #[error("Phase quenching requires a connected vacuum graph; found {components} components")]
    DisconnectedGraph { components: usize },
    #[error("Phase quenching requires a vacuum graph; edge {edge} is not internal")]
    ExternalEdge { edge: EdgeIndex },
    #[error("Phase quenching supports Dirac quarks; edge {edge} is another fermion")]
    UnsupportedFermion { edge: EdgeIndex },
    #[error("Phase quenching does not support a chemical potential on bosonic edge {edge}")]
    BosonicChemicalPotential { edge: EdgeIndex },
    #[error(
        "Quark lines must form unambiguous closed loops; vertex {vertex} has {degree} quark half-edges"
    )]
    QuarkDegree { vertex: usize, degree: usize },
    #[error("Quark flow is not conserved at vertex {vertex}")]
    QuarkFlow { vertex: usize },
    #[error("Quark flavor changes at vertex {vertex}")]
    QuarkFlavor { vertex: usize },
    #[error(transparent)]
    Distribution(#[from] ThermalDistributionCallError),
    #[error("Thermal distribution refers to edge {edge}, absent from the original graph")]
    UnknownEdge { edge: EdgeIndex },
    #[error(
        "Expected a polynomial thermal weight containing only distributions and numbers, got {0}"
    )]
    InvalidWeight(Atom),
    #[error("Expected a thermal_weight call with one argument, got {0}")]
    InvalidWrapper(Atom),
    #[error("Phase quenching requires one complete thermal weight per term, got {0}")]
    IncompleteWeight(Atom),
}

/// Immutable original-graph membership for independent chemical-potential signs.
///
/// Each closed quark loop gets one sign, even when different loops share a
/// flavor. Retain this map while UV contractions remove edges: a surviving
/// distribution still uses its original physical edge ID.
#[derive(Debug)]
pub(crate) struct PhaseQuenchedDifference {
    edge_loops: BTreeMap<EdgeIndex, Option<usize>>,
    loop_count: usize,
}

impl PhaseQuenchedDifference {
    /// Validate a connected vacuum graph and identify its closed Dirac-quark loops.
    pub(crate) fn new(graph: &Graph) -> Result<Self, PhaseQuenchedDifferenceError> {
        let components = graph.count_connected_components(&graph.full_filter());
        if components != 1 {
            return Err(PhaseQuenchedDifferenceError::DisconnectedGraph { components });
        }

        let mut quarks: SuBitGraph = graph.empty_subgraph();
        let mut edge_loops = BTreeMap::new();
        for (pair, edge, data) in graph.underlying.iter_edges() {
            if !matches!(pair, HedgePair::Paired { .. }) {
                return Err(PhaseQuenchedDifferenceError::ExternalEdge { edge });
            }
            edge_loops.insert(edge, None);
            if data.data.is_fermion() {
                let particle = data
                    .data
                    .particle()
                    .expect("a fermion edge has particle metadata");
                if particle.spin != 2
                    || particle.color.abs() != 3
                    || particle.is_self_antiparticle()
                {
                    return Err(PhaseQuenchedDifferenceError::UnsupportedFermion { edge });
                }
                quarks.add(pair);
            } else if data.data.chemical_potential_atom().is_some() {
                return Err(PhaseQuenchedDifferenceError::BosonicChemicalPotential { edge });
            }
        }

        for (vertex, _, _) in graph.underlying.iter_nodes() {
            let incident = graph
                .iter_crown(vertex)
                .filter(|hedge| quarks.includes(hedge))
                .collect::<Vec<_>>();
            if incident.is_empty() {
                continue;
            }
            let [first, second] = incident.as_slice() else {
                return Err(PhaseQuenchedDifferenceError::QuarkDegree {
                    vertex: vertex.into(),
                    degree: incident.len(),
                });
            };
            if graph.superficial_hedge_orientation(*first)
                == graph.superficial_hedge_orientation(*second)
            {
                return Err(PhaseQuenchedDifferenceError::QuarkFlow {
                    vertex: vertex.into(),
                });
            }
            if graph
                .get_edge_data(*first)
                .particle_pdg_code()
                .map(isize::abs)
                != graph
                    .get_edge_data(*second)
                    .particle_pdg_code()
                    .map(isize::abs)
            {
                return Err(PhaseQuenchedDifferenceError::QuarkFlavor {
                    vertex: vertex.into(),
                });
            }
        }

        // Different closed quark cycles carry independent conserved loop
        // circulations. An invertible change of loop basis cannot identify
        // their propagator momenta. Reduced cographs retain original edge IDs
        // and only discard hard coordinates which vanish on every surviving
        // edge, so common-pole reduction cannot merge these sign groups.
        let loops = graph.connected_components(&quarks);
        for (loop_id, quark_loop) in loops.iter().enumerate() {
            for hedge in quark_loop.included_iter() {
                let edge = graph.underlying[&hedge];
                // Literal-zero potentials have no sign dependence. Adjustable
                // named parameters remain active even when their value is zero.
                if graph[edge].chemical_potential_atom().is_some() {
                    edge_loops.insert(edge, Some(loop_id));
                }
            }
        }
        Ok(Self {
            edge_loops,
            loop_count: loops.len(),
        })
    }

    /// Replace each completed weight by `P_PQ(W) - W`, preserving its wrapper.
    ///
    /// This is the `p_PQ - p_QCD` convention. Apply exactly once to each newly
    /// generated CFF: the difference satisfies `D² = -D`, not idempotence.
    /// Expansion is confined to weights; the graph numerator stays factorized.
    ///
    /// Partial hard-edge vacuum replacement commutes with every loop sign
    /// flip: `N(e,d,t,s,o) -> v_d(s)` loses the chemical orientation `o`.
    /// Hence it commutes with this linear difference while original loop labels
    /// are retained. The difference itself is not multiplicative; it must act
    /// on the complete weight, including factors outside the hard subgraph.
    pub(crate) fn apply(&self, expression: &Atom) -> Result<Atom, PhaseQuenchedDifferenceError> {
        self.project_expression(expression.as_view())
    }

    fn project_expression(
        &self,
        expression: AtomView<'_>,
    ) -> Result<Atom, PhaseQuenchedDifferenceError> {
        // This also removes chemical-independent terms in mixed sums and the
        // empty cograph's unit coefficient.
        if !expression.contains_symbol(GS.thermal_weight_wrapper) {
            return Ok(Atom::Zero);
        }
        match expression {
            AtomView::Fun(function) if function.get_symbol() == GS.thermal_weight_wrapper => {
                if function.get_nargs() != 1 {
                    return Err(PhaseQuenchedDifferenceError::InvalidWrapper(
                        expression.to_owned(),
                    ));
                }
                Ok(function!(
                    GS.thermal_weight_wrapper,
                    self.project_weight(function.get(0))?
                ))
            }
            AtomView::Add(sum) => sum.iter().map(|term| self.project_expression(term)).sum(),
            AtomView::Mul(product) => {
                let mut weight_factors = product
                    .iter()
                    .filter(|factor| factor.contains_symbol(GS.thermal_weight_wrapper));
                let weight = weight_factors.next().expect("a thermal wrapper is present");
                if weight_factors.next().is_some() {
                    return Err(PhaseQuenchedDifferenceError::IncompleteWeight(
                        expression.to_owned(),
                    ));
                }
                let projected = self.project_expression(weight)?;
                Ok(product
                    .iter()
                    .filter(|factor| *factor != weight)
                    .fold(projected, |result, factor| result * factor))
            }
            _ => Err(PhaseQuenchedDifferenceError::IncompleteWeight(
                expression.to_owned(),
            )),
        }
    }

    fn validate_weight(&self, weight: AtomView<'_>) -> Result<(), PhaseQuenchedDifferenceError> {
        match weight {
            AtomView::Num(_) => Ok(()),
            AtomView::Fun(function) if function.get_symbol() == GS.thermal_distribution => {
                let call = ThermalDistributionCall::try_from(weight)?;
                if self.edge_loops.contains_key(&call.edge) {
                    Ok(())
                } else {
                    Err(PhaseQuenchedDifferenceError::UnknownEdge { edge: call.edge })
                }
            }
            AtomView::Add(sum) => sum.iter().try_for_each(|term| self.validate_weight(term)),
            AtomView::Mul(product) => product
                .iter()
                .try_for_each(|factor| self.validate_weight(factor)),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                if usize::try_from(exponent).is_err() {
                    return Err(PhaseQuenchedDifferenceError::InvalidWeight(
                        weight.to_owned(),
                    ));
                }
                self.validate_weight(base)
            }
            _ => Err(PhaseQuenchedDifferenceError::InvalidWeight(
                weight.to_owned(),
            )),
        }
    }

    fn project_weight(&self, weight: AtomView<'_>) -> Result<Atom, PhaseQuenchedDifferenceError> {
        self.validate_weight(weight)?;
        // Expose vacuum constants before taking the difference, so they cancel
        // symbolically rather than swallowing exponentially small occupations.
        // N(-1) is the bare signed occupation and derivatives have no vacuum part.
        let original = weight.replace_multiple(&*OCCUPATION_FORM).expand();
        let mut average = original.to_owned();
        for loop_id in 0..self.loop_count {
            let flipped = average.replace_map(|view, _, out| {
                let AtomView::Fun(function) = view else {
                    return;
                };
                if function.get_symbol() != GS.thermal_distribution {
                    return;
                }
                let call = ThermalDistributionCall::try_from(view)
                    .expect("the thermal weight has been validated");
                if self.edge_loops[&call.edge] == Some(loop_id) {
                    **out = ThermalDistributionCall::atom(
                        call.edge.0 as i64,
                        call.derivative_order as i64,
                        call.temperature_flag,
                        call.thermal_sign,
                        -call.orientation_sign,
                    );
                }
            });
            // The independent two-point averages commute and avoid an integer
            // bit mask or a machine-size limit on the number of quark loops.
            average = ((average + flipped) / Atom::num(2)).expand();
        }
        Ok((average - original).expand())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::parse::from_dot::IntoGraph;
    use crate::initialisation;
    use crate::utils;
    use eyre::Result as EyreResult;
    use three_dimensional_reps::ThermalDistributionFactor;

    fn projector() -> PhaseQuenchedDifference {
        PhaseQuenchedDifference {
            edge_loops: BTreeMap::from([
                (EdgeIndex(0), Some(0)),
                (EdgeIndex(1), Some(0)),
                (EdgeIndex(2), Some(1)),
                (EdgeIndex(3), None),
            ]),
            loop_count: 2,
        }
    }

    fn distribution(edge: i64, derivative: i64, contour: i64, orientation: i64) -> Atom {
        ThermalDistributionCall::atom(edge, derivative, 1, contour, orientation)
    }

    fn wrap(weight: Atom) -> Atom {
        function!(GS.thermal_weight_wrapper, weight)
    }

    #[test]
    fn same_loop_factors_share_one_sign() {
        let weight = distribution(0, 0, -1, 1) * distribution(1, 0, -1, -1);
        let reversed = distribution(0, 0, -1, -1) * distribution(1, 0, -1, 1);
        assert_eq!(
            projector().apply(&wrap(weight.to_owned())).unwrap(),
            wrap(((reversed - weight) / Atom::num(2)).expand()),
        );
    }

    #[test]
    fn distinct_loops_are_averaged_independently() {
        let weight = distribution(0, 0, -1, 1) * distribution(2, 0, -1, 1);
        let expected = ((distribution(0, 0, -1, 1) + distribution(0, 0, -1, -1))
            * (distribution(2, 0, -1, 1) + distribution(2, 0, -1, -1))
            / Atom::num(4)
            - &weight)
            .expand();
        assert_eq!(projector().apply(&wrap(weight)).unwrap(), wrap(expected));
    }

    #[test]
    fn derivative_calls_flip_only_the_chemical_orientation() {
        for order in 1..=4 {
            for contour in [-1, 0, 1] {
                let weight = distribution(0, order, contour, 1);
                let expected = (distribution(0, order, contour, -1) - &weight) / Atom::num(2);
                assert_eq!(
                    projector().apply(&wrap(weight)).unwrap(),
                    wrap(expected.expand()),
                );
            }
        }
    }

    #[test]
    fn chemical_independent_weights_cancel_and_numerators_stay_factorized() {
        let projection = projector();
        let boson = distribution(3, 0, 1, 1);
        assert_eq!(
            projection.apply(&wrap(boson.to_owned())).unwrap(),
            Atom::Zero
        );
        assert_eq!(projection.apply(&wrap(Atom::one())).unwrap(), Atom::Zero);
        assert_eq!(projection.apply(&Atom::one()).unwrap(), Atom::Zero);
        let numerator = parse!("(numerator_a + numerator_b)*(numerator_c + numerator_d)");
        let weight = boson * distribution(0, 0, -1, 1);
        assert_eq!(
            projection
                .apply(&(&numerator * wrap(weight.to_owned())))
                .unwrap(),
            numerator * projection.apply(&wrap(weight)).unwrap(),
        );
    }

    #[test]
    fn vacuum_constants_cancel_before_evaluation() {
        let weight = distribution(0, 0, 1, 1);
        let expected = (distribution(0, 0, -1, -1) - distribution(0, 0, -1, 1)) / Atom::num(2);
        assert_eq!(
            projector().apply(&wrap(weight)).unwrap(),
            wrap(expected.expand())
        );
    }

    #[test]
    fn difference_is_not_idempotent() {
        let projection = projector();
        let weight = wrap(distribution(0, 0, -1, 1) * distribution(2, 2, 1, -1));
        let once = projection.apply(&weight).unwrap();
        let twice = projection.apply(&once).unwrap();
        assert_ne!(once, twice);
        let AtomView::Fun(function) = once.as_view() else {
            panic!("nonzero projected weight retains its wrapper");
        };
        assert_eq!(twice, wrap((-function.get(0).to_owned()).expand()));
    }

    #[test]
    fn partial_vacuum_limits_commute_with_correlated_loop_projection() {
        let projection = projector();
        let vacuum_limit = |expression: &Atom, hard_mask: usize| {
            let limited = expression.replace_map(|view, _, out| {
                let AtomView::Fun(function) = view else {
                    return;
                };
                if function.get_symbol() != GS.thermal_distribution {
                    return;
                }
                let call = ThermalDistributionCall::try_from(view).unwrap();
                if hard_mask & (1 << call.edge.0) != 0 {
                    **out = ThermalDistributionFactor::vacuum_atom(
                        call.derivative_order,
                        call.thermal_sign,
                    );
                }
            });
            match limited.as_view() {
                AtomView::Fun(function) if function.get_symbol() == GS.thermal_weight_wrapper => {
                    wrap(function.get(0).expand())
                }
                _ => {
                    assert!(limited.is_zero());
                    limited
                }
            }
        };

        for temperature in [0, 1] {
            for derivative in 0..=3 {
                let factor = |edge, order, contour, orientation| {
                    ThermalDistributionCall::atom(edge, order, temperature, contour, orientation)
                };
                let weight = wrap(
                    (Atom::num(2) + factor(0, 0, 1, 1))
                        * factor(1, 0, 0, -1)
                        * (factor(2, derivative, 1, 1) + factor(2, 0, 1, -1))
                        * (Atom::num(3) + factor(3, derivative, 1, -1))
                        + factor(0, derivative, -1, -1) * factor(2, 0, 0, 1),
                );
                let projected = projection.apply(&weight).unwrap();
                assert!(!projected.is_zero());

                // Every subset includes cuts through loop {0,1}, complete
                // loops, their union, the bosonic edge, and the overall limit.
                for hard_mask in 0..16 {
                    let early = vacuum_limit(&projected, hard_mask);
                    let late = projection.apply(&vacuum_limit(&weight, hard_mask)).unwrap();
                    assert_eq!(
                        early, late,
                        "temperature={temperature}, derivative={derivative}, hard_mask={hard_mask}"
                    );
                    if hard_mask & 0b0111 == 0b0111 {
                        assert!(
                            early.is_zero(),
                            "all chemical dependence is in the vacuum limit"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn mixed_sums_drop_vacuum_terms_and_reject_separate_weight_products() {
        let projection = projector();
        let weight = wrap(distribution(0, 0, -1, 1));
        let coefficient = parse!("coefficient_a + coefficient_b");
        assert_eq!(
            projection
                .apply(&(&coefficient * (&weight + Atom::one())))
                .unwrap(),
            coefficient * projection.apply(&weight).unwrap(),
        );
        assert!(matches!(
            projection.apply(&(weight * wrap(distribution(1, 0, -1, 1)))),
            Err(PhaseQuenchedDifferenceError::IncompleteWeight(_)),
        ));
    }

    #[test]
    fn rejects_unknown_edges_and_nonpolynomial_weights() {
        let projection = projector();
        assert!(matches!(
            projection.apply(&wrap(distribution(9, 0, 1, 1))),
            Err(PhaseQuenchedDifferenceError::UnknownEdge { .. }),
        ));
        assert!(matches!(
            projection.apply(&wrap(distribution(0, 0, 1, 1).pow(-1))),
            Err(PhaseQuenchedDifferenceError::InvalidWeight(_)),
        ));
    }

    #[test]
    fn original_graph_distinguishes_same_flavor_loops() -> EyreResult<()> {
        initialisation::test_initialise()?;
        let model = utils::load_generic_model("sm");
        let graph: Graph = r#"
            digraph two_quark_loops {
                node [num=1]; edge [num=1];
                A -> B [id=0 particle="d"];
                B -> A [id=1 particle="d"];
                C -> D [id=2 particle="d"];
                D -> C [id=3 particle="d"];
                A -> C [id=4 particle="g"];
            }
        "#
        .into_graph(&model)?;
        let projection = PhaseQuenchedDifference::new(&graph)?;
        assert_eq!(projection.loop_count, 2);
        assert_eq!(
            projection.edge_loops[&EdgeIndex(0)],
            projection.edge_loops[&EdgeIndex(1)]
        );
        assert_eq!(
            projection.edge_loops[&EdgeIndex(2)],
            projection.edge_loops[&EdgeIndex(3)]
        );
        assert_ne!(
            projection.edge_loops[&EdgeIndex(0)],
            projection.edge_loops[&EdgeIndex(2)]
        );
        assert_eq!(projection.edge_loops[&EdgeIndex(4)], None);
        Ok(())
    }

    #[test]
    fn rejects_ambiguous_quark_intersections() -> EyreResult<()> {
        initialisation::test_initialise()?;
        let model = utils::load_generic_model("sm");
        let graph: Graph = r#"
            digraph ambiguous_loops {
                node [num=1]; edge [num=1];
                A -> A [id=0 particle="d"];
                A -> A [id=1 particle="d"];
            }
        "#
        .into_graph(&model)?;
        assert!(matches!(
            PhaseQuenchedDifference::new(&graph),
            Err(PhaseQuenchedDifferenceError::QuarkDegree { degree: 4, .. }),
        ));
        Ok(())
    }

    #[test]
    fn pure_bosonic_graph_has_zero_difference() -> EyreResult<()> {
        initialisation::test_initialise()?;
        let model = utils::load_generic_model("sm");
        let graph: Graph = r#"
            digraph gluon_loop {
                node [num=1]; edge [num=1];
                A -> A [id=0 particle="g"];
            }
        "#
        .into_graph(&model)?;
        let projection = PhaseQuenchedDifference::new(&graph)?;
        assert_eq!(projection.loop_count, 0);
        assert_eq!(
            projection.apply(&wrap(distribution(0, 0, 1, 1)))?,
            Atom::Zero
        );
        Ok(())
    }

    #[test]
    fn rejects_open_quark_chains_flavor_changes_and_flow_mismatches() -> EyreResult<()> {
        initialisation::test_initialise()?;
        let model = utils::load_generic_model("sm");
        let graphs = [
            (
                r#"digraph open_chain { node [num=1]; edge [num=1];
                    A -> B [id=0 particle="d"]; B -> A [id=1 particle="g"]; }"#,
                "quark half-edges",
            ),
            (
                r#"digraph flavor_change { node [num=1]; edge [num=1];
                    A -> B [id=0 particle="d"]; B -> A [id=1 particle="u"]; }"#,
                "flavor changes",
            ),
            (
                r#"digraph flow_mismatch { node [num=1]; edge [num=1];
                    A -> B [id=0 particle="d"]; A -> B [id=1 particle="d"]; }"#,
                "flow is not conserved",
            ),
        ];
        for (source, message) in graphs {
            let graph: Graph = source.into_graph(&model)?;
            let error = PhaseQuenchedDifference::new(&graph).unwrap_err();
            assert!(error.to_string().contains(message), "{error}");
        }
        Ok(())
    }
}
