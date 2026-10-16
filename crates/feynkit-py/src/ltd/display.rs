//! Inspection data retains native coefficients, multiplicities and affine maps.
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Write;

use feynkit_cff::generalized::{HybridSurfaceID, LinearEnergyExpr, LinearSurfaceKind};
use linnet::half_edge::involution::EdgeIndex;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use serde_json::{Value, json};
use symbolica::{
    api::python::PythonExpression,
    atom::{Atom, AtomCore, Symbol},
    printer::PrintOptions,
};

use super::PyLtdRepresentation;
use crate::display::{energy_graph_svg, escape_html};

impl PyLtdRepresentation {
    fn graph_labels(&self) -> PyResult<String> {
        let mut labels = Vec::new();
        for edge in self.indices.internal.values() {
            let momentum = feynkit_graph::symbols::momentum()
                .call(*edge)
                .printer(PrintOptions::typst())
                .to_string();
            let component = feynkit_cff::symbols::external_energy_atom(EdgeIndex(*edge))
                .printer(PrintOptions::typst())
                .to_string();
            let energy = feynkit_cff::symbols::on_shell_atom(EdgeIndex(*edge))
                .printer(PrintOptions::typst())
                .to_string();
            labels.push((format!("q-{edge}"), format!("[$ {momentum} $]")));
            for (name, sign) in [("plus", "+"), ("minus", "-")] {
                labels.push((
                    format!("q-{edge}-{name}"),
                    format!("[$ {momentum} thin {sign} $]"),
                ));
                labels.push((
                    format!("pole-{edge}-{name}"),
                    format!("align(center, stack(spacing: 2pt, box[$ {component} $], box[$ = {sign} {energy} $]))"),
                ));
            }
        }
        for edge in &self.parsed.external_edges {
            let coordinate = edge.edge_id - self.indices.internal.len();
            let physical = self.indices.external[&coordinate];
            let momentum = feynkit_graph::symbols::external_momentum()
                .call(coordinate)
                .printer(PrintOptions::typst())
                .to_string();
            labels.push((format!("p-{physical}"), format!("[$ {momentum} $]")));
        }
        // Compile all states once with Linnet's label compiler and glyph merger.
        // Hover and pair selection then only change a local SVG reference.
        let mut svg = String::from("<svg data-ltd-labels><defs>");
        if !labels.is_empty() {
            let source = format!(
                "#set page(width: auto, height: auto, margin: 0pt, fill: none)\n\
                 #set text(size: 9pt)\n{}",
                labels
                    .iter()
                    .map(|(_, source)| format!("#{source}"))
                    .collect::<Vec<_>>()
                    .join("\n#pagebreak()\n"),
            );
            let pages = typst_renderer::Document::compile_sources(
                &BTreeMap::from([("main.typ".into(), source.into_bytes())]),
                "svg",
            )
            .map_err(PyRuntimeError::new_err)?
            .into_iter()
            .map(|page| String::from_utf8(page).map_err(PyRuntimeError::new_err))
            .collect::<PyResult<Vec<_>>>()?;
            let typeset =
                linnest::svg::Typeset::read(&pages, false).map_err(PyRuntimeError::new_err)?;
            for ((key, _), page) in labels.iter().zip(&typeset.pages) {
                write!(svg, "<g id=\"ltd-label-{key}\" data-label-key=\"{key}\" data-width=\"{}\" data-height=\"{}\">{}</g>", page.width, page.height, page.body).unwrap();
            }
            svg.push_str(&typeset.defs);
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

    pub(super) fn explorer_html(&self, py: Python<'_>, residue: Option<usize>) -> PyResult<String> {
        if self.drawing.get().is_none() {
            let mut drawing = energy_graph_svg(
                py,
                &self.diagram,
                &feynkit_graph::SceneOptions {
                    momentum_arrows: true,
                    // Leave room for the pole chevrons on the propagator itself.
                    momentum_arrow_offset: 0.7,
                    ..Default::default()
                },
            )?;
            drawing.push_str(&self.graph_labels()?);
            let _ = self.drawing.set(drawing);
        }
        let data = self.explorer_data(py, residue)?;
        let data = data.to_string().replace('<', "\\u003c");
        Ok(format!(
            "<style>{}{}{}</style><section class=\"feynkit-ltd-result feynkit-expression\" data-feynkit-notebook data-linnet-frame-owner aria-label=\"{} LTD\"><div class=\"ltd-product\"></div><template data-drawing>{}</template><script type=\"application/json\" data-ltd>{data}</script><noscript>Enable notebook JavaScript to inspect LTD residues, or use to_expression().</noscript></section><script>{}</script><script>{}</script><script>{}</script>",
            spynso3::display::NOTEBOOK_STYLE,
            include_str!("../expression.css"),
            include_str!("display.css"),
            escape_html(self.diagram.name()),
            self.drawing.get().expect("rendered above"),
            include_str!("../expression.js"),
            include_str!("display.js"),
            include_str!("../notebook.js"),
        ))
    }
}
