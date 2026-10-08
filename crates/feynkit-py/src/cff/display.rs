//! CFF presentation owns only inspection state; native trees and graph IDs own the data.
use std::collections::BTreeMap;

use pyo3::prelude::*;
use serde_json::{Value, json};
use symbolica::{
    api::python::PythonExpression,
    atom::{Atom, AtomCore},
};

use super::PyCffRepresentation;
use crate::display::{energy_graph_svg, escape_html};

#[derive(Clone, Copy)]
pub(super) enum CffScope {
    Representation,
    Orientation(usize),
    Family { orientation: usize, family: usize },
}

impl PyCffRepresentation {
    fn explorer_data(&self, py: Python<'_>, scope: CffScope) -> PyResult<Value> {
        // Energy and H indices occupy separate arenas. Give the browser one dense
        // key, while keeping the native category and index for the displayed label.
        let surfaces = self.surfaces();
        let selected_orientation = match scope {
            CffScope::Representation => None,
            CffScope::Orientation(id)
            | CffScope::Family {
                orientation: id, ..
            } => Some(id),
        };
        let selected_family = match scope {
            CffScope::Family { family, .. } => Some(family),
            _ => None,
        };
        let mut math_cache = BTreeMap::<Atom, String>::new();
        let format_text = py
            .import("symbolica.community.tensor")?
            .getattr("format_tensor")?;
        let mut math = |expr: Atom| -> PyResult<String> {
            if expr == Atom::num(1) {
                return Ok("1".into());
            }
            if let Some(html) = math_cache.get(&expr) {
                return Ok(html.clone());
            }
            let html =
                crate::display::expression_html(py, PythonExpression { expr: expr.clone() })?;
            math_cache.insert(expr, html.clone());
            Ok(html)
        };
        let mut orientations = Vec::new();
        for orientation in self
            .orientations()
            .into_iter()
            .filter(|o| selected_orientation.is_none_or(|id| o.id() == id))
        {
            let mut family_ids = Vec::new();
            let mut terms = Vec::new();
            let mut contributions = Vec::new();
            for (id, term) in orientation
                .inner
                .terms
                .iter()
                .enumerate()
                .filter(|(id, _)| selected_family.is_none_or(|selected| selected == *id))
            {
                family_ids.push(id);
                terms.push(term.path.clone());
                let numerator = term.numerator_atom(&surfaces);
                contributions.push(json!({"coefficient": term.coefficient.to_canonical_string(),
                    "energies": term.energies, "numerator": numerator.to_canonical_string(),
                    "numerator_html": math(numerator)?,
                    "energy_map": term.energy_map.iter().map(|(e,a)| Ok((*e,format_text.call1((PythonExpression { expr: a.clone() },))?.extract::<String>()?))).collect::<PyResult<BTreeMap<_,_>>>()?,
                    "origin": term.origin,
                }));
            }
            orientations.push(json!({"id": orientation.id(), "directions": orientation.edge_orientations().into_iter().collect::<BTreeMap<_, _>>(),
                "terms": terms, "family_ids": family_ids, "contributions": contributions}));
        }
        let surface_math = surfaces
            .iter()
            .map(|s| Ok((s.atom(false), math(s.atom(true))?)))
            .collect::<PyResult<BTreeMap<_, _>>>()?;
        let mut energy_html = self
            .diagram
            .edges()
            .filter(|(_, _, edge)| edge.external.is_none())
            .map(|(edge, _, _)| {
                Ok((
                    edge.0.to_string(),
                    math(feynkit_cff::symbols::on_shell_atom(
                        linnet::half_edge::involution::EdgeIndex(edge.0),
                    ))?,
                ))
            })
            .collect::<PyResult<BTreeMap<_, _>>>()?;
        energy_html.insert(
            "e".into(),
            math(feynkit_cff::symbols::on_shell().call(symbolica::symbol!("feynkit_display::e")))?,
        );
        Ok(json!({
            "generalized": self.generalized,
            "normalization": self.normalization.to_canonical_string(),
            "normalization_html": if self.generalized { math(self.normalization.clone())? } else { String::new() },
            "name": self.diagram.name(),
            "loops": self.diagram.loop_count(),
            "energy_edges": self.energy_edges,
            "energy_html": energy_html,
            "scope": match scope { CffScope::Representation => "representation", CffScope::Orientation(_) => "orientation", CffScope::Family {..} => "family" },
            "pole_order": self.pole_order,
            "orientations": orientations,
            "surfaces": surfaces.iter().map(|s| json!({
                "index": s.index(),
                "origin": s.origin(), "numerator_only": s.numerator_only(),
                "support": s.energy_coefficients().iter().filter_map(|(e,c)| (!c.expr.is_zero()).then_some(*e)).collect::<Vec<_>>(),
                "symbol": s.atom(false).to_canonical_string(),
                "expression": s.atom(true).to_canonical_string(),
                "expression_html": surface_math.get(&s.atom(false)),
                "kind": if s.kind() == "H" { "h" } else { "energy" },
                "e": s.energy_coefficients().iter().filter_map(|(e,c)| (c.expr == Atom::num(1)).then_some(*e)).collect::<Vec<_>>(),
                "negative": s.energy_coefficients().iter().filter_map(|(e,c)| (c.expr == Atom::num(-1)).then_some(*e)).collect::<Vec<_>>(),
                "q": s.external_shift().iter().filter_map(|(e,c)| i64::try_from(c.expr.as_view()).ok().map(|c| (*e,c))).collect::<Vec<_>>(), "v": s.vertices(),
            })).collect::<Vec<_>>(),
            "edges": self.diagram.edges().map(|(id, ends, edge)| json!({
                "id": id.0, "source": ends.source.map(|v| v.0), "target": ends.target.map(|v| v.0),
                "external": edge.external.is_some(),
            })).collect::<Vec<_>>(),
        }))
    }

    pub(super) fn explorer_html(&self, py: Python<'_>, scope: CffScope) -> PyResult<String> {
        if self.drawing.get().is_none() {
            let svg = energy_graph_svg(py, &self.diagram, &Default::default())?;
            let _ = self.drawing.set(svg);
        }
        // JSON is inert data, but a literal closing script tag still ends an HTML
        // raw-text element. Escape '<' even in graph names supplied by the user.
        let data = self.explorer_data(py, scope)?;
        let kind = match scope {
            CffScope::Representation => "CffRepresentation",
            CffScope::Orientation(_) => "CffOrientation",
            CffScope::Family { .. } => "CrossFreeFamily",
        };
        let data = data.to_string().replace('<', "\\u003c");
        Ok(format!(
            "<style>{}{}{}</style><section class=\"feynkit-cff-result feynkit-expression\" data-feynkit-notebook data-linnet-frame-owner><header class=\"hs-meta\"><code title=\"{}\">{kind}</code><span data-metadata></span></header>{}<template data-drawing>{}</template><script type=\"application/json\" data-cff>{data}</script><noscript>Enable notebook JavaScript to explore orientations and families. Use to_expression() for the symbolic expression.</noscript></section><script>{}</script><script>{}</script><script>{}</script>",
            spynso3::display::NOTEBOOK_STYLE,
            include_str!("../expression.css"),
            include_str!("display.css"),
            escape_html(self.diagram.name()),
            include_str!("display.html"),
            self.drawing.get().expect("rendered above"),
            include_str!("../expression.js"),
            include_str!("display.js"),
            include_str!("../notebook.js"),
        ))
    }
}
