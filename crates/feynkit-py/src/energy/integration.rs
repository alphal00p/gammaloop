//! Shared physical-edge adapter for native three-dimensional generation.
use crate::error;
use feynkit_cff::generalized::{
    LinearEnergyExpr, MomentumSignature, ParsedGraph,
    graph_io::{EnergyEdgeIndexMap, ParsedGraphExternalEdge, ParsedGraphInternalEdge},
};
use feynkit_graph::FeynmanDiagram;
use pyo3::{exceptions::PyValueError, prelude::*};
use std::{collections::BTreeMap, sync::Arc};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::atom::AtomField,
    poly::PolyVariable,
};

pub(crate) struct EnergyGraph {
    pub parsed: ParsedGraph,
    pub indices: EnergyEdgeIndexMap,
}
impl EnergyGraph {
    pub fn new(diagram: &FeynmanDiagram) -> PyResult<Self> {
        let basis = diagram.loop_momentum_basis();
        let internal_count = diagram
            .edges()
            .filter(|(_, _, e)| e.external.is_none())
            .count();
        let mut internal_edges = Vec::new();
        let mut external_edges = Vec::new();
        let mut internal = BTreeMap::new();
        for (id, ends, edge) in diagram.edges() {
            if edge.is_dummy {
                return Err(error::DiagramError::new_err(
                    "energy integration requires an ordinary diagram without auxiliary cut edges",
                ));
            }
            let signature = &basis.edge_signatures[&id];
            let external_signature = signature
                .external
                .integer_coefficients()
                .into_iter()
                .map(|c| c as i32)
                .collect();
            if edge.external.is_some() {
                let coordinate = basis
                    .external_edges
                    .iter()
                    .position(|e| *e == id)
                    .ok_or_else(|| {
                        error::DiagramError::new_err(
                            "external edge is absent from the momentum basis",
                        )
                    })?;
                external_edges.push(ParsedGraphExternalEdge {
                    edge_id: internal_count + coordinate,
                    source: ends.source.map(|v| v.0),
                    destination: ends.target.map(|v| v.0),
                    label: format!("p{coordinate}"),
                    external_coefficients: external_signature,
                });
            } else {
                let local = internal_edges.len();
                internal.insert(local, id.0);
                internal_edges.push(ParsedGraphInternalEdge {
                    edge_id: local,
                    tail: ends
                        .source
                        .ok_or_else(|| error::DiagramError::new_err("internal edge has no source"))?
                        .0,
                    head: ends
                        .target
                        .ok_or_else(|| error::DiagramError::new_err("internal edge has no target"))?
                        .0,
                    label: format!("q{}", id.0),
                    mass_key: Some(
                        diagram
                            .model()
                            .particle_by_id(edge.particle)
                            .map_err(error::model)?
                            .symbolic_mass(diagram.model())
                            .to_canonical_string(),
                    ),
                    signature: MomentumSignature {
                        loop_signature: signature
                            .loops
                            .integer_coefficients()
                            .into_iter()
                            .map(|c| c as i32)
                            .collect(),
                        external_signature,
                    },
                    had_pow: false,
                });
            }
        }
        let parsed = ParsedGraph {
            internal_edges,
            external_edges,
            initial_state_cut_edges: Vec::new(),
            loop_names: (0..basis.loop_edges.len())
                .map(|i| format!("k{i}"))
                .collect(),
            external_names: (0..basis.external_edges.len())
                .map(|i| format!("p{i}"))
                .collect(),
            node_name_to_internal: diagram
                .vertices()
                .map(|(id, _)| (format!("v{}", id.0), id.0))
                .collect(),
        };
        let indices = EnergyEdgeIndexMap {
            internal,
            external: basis
                .external_edges
                .iter()
                .enumerate()
                .map(|(i, id)| (i, id.0))
                .collect(),
            orientation_edge_count: diagram
                .edges()
                .map(|(id, _, _)| id.0 + 1)
                .max()
                .unwrap_or(0),
        };
        Ok(Self { parsed, indices })
    }

    pub fn physical_energy(&self, expression: &LinearEnergyExpr) -> LinearEnergyExpr {
        expression
            .clone()
            .remap_energy_edges(&self.indices.internal, &self.indices.external)
    }

    pub fn numerator_at(&self, numerator: &Atom, map: &[LinearEnergyExpr]) -> Atom {
        numerator.replace_multiple(self.indices.internal.iter().map(|(local, physical)| {
            symbolica::id::Replacement::new(
                feynkit_cff::symbols::external_energy_atom(
                    linnet::half_edge::involution::EdgeIndex(*physical),
                ),
                self.physical_energy(&map[*local]).to_atom(&[]),
            )
        }))
    }

    // Treat edge energies as independent polynomial coordinates: routing before
    // bounding would erase the edge-local degree certificate used by the generator.
    pub fn numerator_bounds(&self, numerator: &Atom) -> PyResult<Vec<(usize, usize)>> {
        for coordinate in 0..self.parsed.loop_names.len() {
            let energy = feynkit_graph::symbols::loop_momentum()
                .call((coordinate, symbolica::symbol!("spenso::cind").call(0)));
            if numerator.contains(energy.as_view()) {
                return Err(PyValueError::new_err(
                    "express numerator energy dependence in physical Q(edge, cind(0)) coordinates before integration",
                ));
            }
        }
        let parameters: Vec<_> = self
            .indices
            .internal
            .values()
            .map(|edge| {
                feynkit_cff::symbols::external_energy_atom(
                    linnet::half_edge::involution::EdgeIndex(*edge),
                )
            })
            .collect();
        if parameters.is_empty() {
            return Ok(Vec::new());
        }
        let requested = parameters
            .iter()
            .cloned()
            .map(PolyVariable::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map_err(PyValueError::new_err)?;
        let field = AtomField {
            statistical_zero_test: false,
            cancel_check_on_division: false,
            custom_normalization: None,
        };
        let polynomial = numerator
            .try_to_polynomial::<_, u32>(&field, Some(Arc::new(requested.clone())))
            .map_err(|_| {
                PyValueError::new_err("numerator must be polynomial in Q(edge, cind(0))")
            })?;
        let indices = requested
            .iter()
            .map(|v| {
                polynomial
                    .variables()
                    .iter()
                    .position(|p| p == v)
                    .ok_or_else(|| {
                        PyValueError::new_err("numerator polynomial changed its energy coordinates")
                    })
            })
            .collect::<PyResult<Vec<_>>>()?;
        let mut bounds = vec![0; parameters.len()];
        for (powers, coefficient) in polynomial.to_multivariate_polynomial_list(&indices, true) {
            let coefficient = coefficient.flatten(false);
            if parameters.iter().any(|p| coefficient.contains(p.as_view())) {
                return Err(PyValueError::new_err(
                    "numerator must be polynomial in Q(edge, cind(0)); energy-dependent functions and denominators are unsupported",
                ));
            }
            for (local, index) in indices.iter().enumerate() {
                bounds[local] = bounds[local].max(powers[*index] as usize);
            }
        }
        Ok(bounds
            .into_iter()
            .enumerate()
            .filter(|(_, degree)| *degree != 0)
            .collect())
    }
}
