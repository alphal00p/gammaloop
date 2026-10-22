//! Inspection data retains native coefficients, multiplicities and affine maps.
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Write;

use feynkit_cff::generalized::{
    HybridSurfaceID, LinearEnergyExpr, LinearSurfaceKind,
    graph_io::{EnergyEdgeIndexMap, ParsedGraphExternalEdge},
};
use linnet::half_edge::involution::EdgeIndex;
use linnet_render_py::PyDiagramView;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use serde_json::{Value, json};
use symbolica::{
    api::python::PythonExpression,
    atom::{Atom, Symbol},
    printer::{PrintMode, PrintOptions, PrintState},
};

use super::PyLtdRepresentation;
use crate::display::{energy_graph_view, escape_html};

impl PyLtdRepresentation {
    fn graph_labels(
        indices: &EnergyEdgeIndexMap,
        external_edges: &[ParsedGraphExternalEdge],
    ) -> PyResult<String> {
        let options = PrintOptions {
            mode: PrintMode::MathML,
            ..PrintOptions::default()
        };
        let mathml = |atom: Atom| {
            spenso::network::tags::mathml::render(
                atom.as_view(),
                &options,
                &PrintState::new(),
                &|_, _| None,
            )
            .ok_or_else(|| PyRuntimeError::new_err("LTD label exceeds native MathML limits"))
        };
        let mut labels = Vec::new();
        for edge in indices.internal.values() {
            let momentum = mathml(feynkit_graph::symbols::momentum().call(*edge))?;
            let component = mathml(feynkit_cff::symbols::external_energy_atom(EdgeIndex(*edge)))?;
            let energy = mathml(feynkit_cff::symbols::on_shell_atom(EdgeIndex(*edge)))?;
            labels.push((format!("q-{edge}"), momentum.clone()));
            for (name, sign) in [("plus", "+"), ("minus", "−")] {
                labels.push((
                    format!("q-{edge}-{name}"),
                    format!("<mrow>{momentum}<mspace width=\"0.167em\"/><mo>{sign}</mo></mrow>"),
                ));
                labels.push((
                    format!("pole-{edge}-{name}"),
                    format!("<mtable rowspacing=\"2px\" columnalign=\"center\"><mtr><mtd>{component}</mtd></mtr><mtr><mtd><mrow><mo>=</mo><mo>{sign}</mo>{energy}</mrow></mtd></mtr></mtable>"),
                ));
            }
        }
        for edge in external_edges {
            let coordinate = edge.edge_id - indices.internal.len();
            let physical = indices.external[&coordinate];
            let momentum = mathml(feynkit_graph::symbols::external_momentum().call(coordinate))?;
            labels.push((format!("p-{physical}"), momentum));
        }
        // Inert fragments share Symbolica/Spenso notation. The browser measures
        // the selected fragment in its visible label slot, including pole rows.
        let mut svg = String::from("<svg data-ltd-labels><defs>");
        for (key, fragment) in labels {
            write!(svg, "<g id=\"ltd-label-{key}\" data-label-key=\"{key}\"><foreignObject width=\"1\" height=\"1\" overflow=\"visible\" pointer-events=\"none\"><div xmlns=\"http://www.w3.org/1999/xhtml\" style=\"display:flex;font-size:9px;line-height:1;width:max-content;white-space:nowrap\"><math xmlns=\"http://www.w3.org/1998/Math/MathML\">{fragment}</math></div></foreignObject></g>").unwrap();
        }
        svg.push_str("</defs></svg>");
        Ok(svg)
    }

    fn linear_data(expression: &LinearEnergyExpr) -> Value {
        json!({
            "internal": expression.internal_terms.iter().map(|(i,c)| (i.0,c.to_string())).collect::<Vec<_>>(),
            "external": expression.external_terms.iter().map(|(i,c)| (i.0,c.to_string())).collect::<Vec<_>>(),
            "constant": expression.constant.to_string(),
        })
    }

    fn explorer_data(&self, py: Python<'_>, residue: Option<usize>) -> PyResult<Value> {
        let cache = &self.inner.expression.surfaces.linear_surface_cache;
        let surface_ids: HashMap<_, _> = cache
            .iter_enumerated()
            .map(|(id, s)| (s.expression.clone(), id.0))
            .collect();
        let repeated = feynkit_cff::generalized::graph_io::repeated_groups(&self.parsed);
        let residues: Vec<_> = self.inner.expression.orientations.iter().enumerate().filter(|(id, _)| residue.is_none_or(|selected| selected == *id)).map(|(id, r)| {
            // Raised poles repeat half-energy factors. Preserve these powers in
            // the expression, but count each actual cut only once in the graph.
            let cuts: BTreeSet<_> = r.variants.iter().flat_map(|v| v.half_edges.iter().map(|e| e.0)).collect();
            let pairs: BTreeMap<_, _> = r.edge_energy_map.iter().enumerate()
                .filter(|(edge,_)| !cuts.contains(edge))
                .filter_map(|(edge, map)| {
                    let minus = (map.clone() + LinearEnergyExpr::ose(EdgeIndex(edge), -1)).canonical();
                    let plus = (map.clone() + LinearEnergyExpr::ose(EdgeIndex(edge), 1)).canonical();
                    Some((edge, [*surface_ids.get(&minus)?, *surface_ids.get(&plus)?]))
                }).collect();
            let terms: Vec<_> = r.variants.iter().map(|v| {
                let chains: Vec<Vec<_>> = v.denominator.get_bottom_layer().into_iter().filter_map(|leaf| {
                    let mut node = v.denominator.get_node(leaf);
                    let mut chain = Vec::new();
                    loop {
                        match node.data {
                            HybridSurfaceID::Linear(id) => chain.push(id.0),
                            HybridSurfaceID::Infinite => return None,
                            HybridSurfaceID::Unit => {},
                            _ => unreachable!("LTD generation interns affine surfaces"),
                        }
                        let Some(parent) = node.parent else { break };
                        node = v.denominator.get_node(parent);
                    }
                    chain.reverse();
                    Some(chain)
                }).collect();
                json!({"coefficient": v.prefactor.to_string(), "energies": v.half_edges.iter().map(|e|e.0).collect::<Vec<_>>(), "chains": chains})
            }).collect();
            json!({
                "id": id, "cuts": cuts, "signs": r.data.orientation,
                "pole_orders": cuts.iter().map(|edge| (*edge, repeated.iter().find(|group| group.edge_ids.contains(edge)).map_or(1, |group| group.edge_ids.len()))).collect::<BTreeMap<_,_>>(),
                "loop_map": r.loop_energy_map.iter().map(Self::linear_data).collect::<Vec<_>>(),
                "edge_map": r.edge_energy_map.iter().map(Self::linear_data).collect::<Vec<_>>(),
                "propagator_factors": pairs, "terms": terms,
            })
        }).collect();
        let edges: Vec<_> = self.indices.internal.values().copied().collect();
        let basis = self.diagram.loop_momentum_basis();
        let energy_html = edges
            .iter()
            .map(|edge| {
                Ok((
                    *edge,
                    crate::display::expression_html(
                        py,
                        PythonExpression {
                            expr: feynkit_cff::symbols::on_shell_atom(EdgeIndex(*edge)),
                        },
                    )?,
                ))
            })
            .collect::<PyResult<BTreeMap<_, _>>>()?;
        let component = symbolica::symbol!("spenso::cind").call(0);
        let momentum_html = |symbol: Symbol, id: usize| {
            crate::display::expression_html(
                py,
                PythonExpression {
                    expr: symbol.call_args([Atom::num(id), component.clone()]),
                },
            )
        };
        let edge_momenta = edges
            .iter()
            .map(|edge| {
                Ok((
                    *edge,
                    momentum_html(feynkit_graph::symbols::momentum(), *edge)?,
                ))
            })
            .collect::<PyResult<BTreeMap<_, _>>>()?;
        let loop_momenta = (0..self.parsed.loop_names.len())
            .map(|id| momentum_html(feynkit_graph::symbols::loop_momentum(), id))
            .collect::<PyResult<Vec<_>>>()?;
        let external_momenta = (0..self.indices.external.len())
            .map(|id| momentum_html(feynkit_graph::symbols::external_momentum(), id))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(json!({
            "name": self.diagram.name(), "loops": self.parsed.loop_names.len(),
            "scope": if residue.is_some() { "residue" } else { "representation" },
            "edges": edges, "routing": self.parsed.internal_edges,
            "energy_html": energy_html,
            "momentum_html": edge_momenta,
            "loop_energy_html": loop_momenta,
            "external_energy_html": external_momenta,
            "external_routing": self.parsed.external_edges.iter().map(|e| (self.indices.external[&(e.edge_id - edges.len())],e)).collect::<BTreeMap<_,_>>(),
            "reference_chords": basis.loop_edges.iter().map(|id| edges.iter().position(|e| *e == id.0)).collect::<Vec<_>>(),
            "residues": residues,
            "surfaces": cache.iter_enumerated().map(|(id,s)| json!({
                "id": id.0, "kind": if s.kind == LinearSurfaceKind::Hsurface { "H" } else { "E" },
                "expression": Self::linear_data(&s.expression),
            })).collect::<Vec<_>>(),
        }))
    }

    pub(crate) fn explorer_view(
        &self,
        py: Python<'_>,
        residue: Option<usize>,
    ) -> PyResult<PyDiagramView> {
        let drawing = energy_graph_view(
            py,
            &self.diagram,
            &feynkit_graph::SceneOptions {
                momentum_arrows: true,
                // Leave room for the pole chevrons on the propagator itself.
                momentum_arrow_offset: 0.7,
                ..Default::default()
            },
        )?;
        let labels = Self::graph_labels(&self.indices, &self.parsed.external_edges)?;
        let data = self.explorer_data(py, residue)?;
        let data = data.to_string().replace('<', "\\u003c");
        Ok(drawing.with_html(format!(
            "<style>{}{}{}</style><section class=\"feynkit-ltd-result feynkit-expression\" data-feynkit-notebook data-linnet-frame-owner aria-label=\"{} LTD\"><div class=\"ltd-product\"></div><template data-drawing><!--linnet-svg-0-->{labels}</template><script type=\"application/json\" data-ltd>{data}</script><noscript>Enable notebook JavaScript to inspect LTD residues, or use to_expression().</noscript></section><script>{}</script><script>{}</script><script>{}</script>",
            spynso3::display::NOTEBOOK_STYLE,
            include_str!("../expression.css"),
            include_str!("display.css"),
            escape_html(self.diagram.name()),
            include_str!("../expression.js"),
            include_str!("display.js"),
            include_str!("../notebook.js"),
        )))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::{EnergyEdgeIndexMap, ParsedGraphExternalEdge, PyLtdRepresentation};

    #[test]
    fn native_label_states_preserve_physical_ids_and_two_row_poles() {
        let indices = EnergyEdgeIndexMap {
            internal: BTreeMap::from([(0, 4)]),
            external: BTreeMap::from([(0, 9)]),
            orientation_edge_count: 10,
        };
        let external = ParsedGraphExternalEdge {
            edge_id: 1,
            source: None,
            destination: Some(0),
            label: "p0".into(),
            external_coefficients: vec![1],
        };
        let labels = PyLtdRepresentation::graph_labels(&indices, &[external]).unwrap();
        assert_eq!(labels.matches("data-label-key=").count(), 6);
        assert_eq!(labels.matches("<mtable ").count(), 2);
        assert_eq!(labels.matches("<mtr>").count(), 4);
        for key in [
            "q-4",
            "q-4-plus",
            "q-4-minus",
            "pole-4-plus",
            "pole-4-minus",
            "p-9",
        ] {
            assert!(labels.contains(&format!("data-label-key=\"{key}\"")));
        }
        assert!(labels.contains("<mo>=</mo><mo>+</mo>"));
        assert!(labels.contains("<mo>=</mo><mo>−</mo>"));
        assert!(labels.contains("<mi>E</mi>"));
        assert!(!labels.contains("data-symbolica-"));
        assert!(!labels.contains("<path"));
        assert!(!labels.contains("<use"));
        // The external key is physical e9, while its displayed momentum uses
        // external coordinate zero, independently of internal edge numbering.
        let external = labels.split("data-label-key=\"p-9\"").nth(1).unwrap();
        assert!(external.contains("<mn>0</mn>"));
        assert!(!external.contains("<mn>9</mn>"));
    }
}
