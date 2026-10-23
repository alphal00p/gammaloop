//! Compact metadata inspectors. Native radio inputs keep slot selection usable
//! in notebook HTML without scripts, frames, or a widget runtime.
use std::sync::atomic::{AtomicU64, Ordering};

use spenso::structure::TensorStructure;
use spenso::structure::dimension::Dimension;
use spenso::structure::partial::PartialIndex;
use spenso::structure::representation::{LibraryRep, Representation};
use spenso::structure::slot::{Slot, SlotMatch, SlotMatcher};

use crate::display;
use crate::display::{DisplaySettings, IndexAliases, TensorDisplayMode, paging};
use crate::metadata::{SpensoRepresentationName, SpensoTensorStructure};
use idenso::tensor::composition;
use pyo3::{PyResult, Python};
use spenso::structure::abstract_index::AIND_SYMBOLS;
use spenso::structure::partial::PartialSlot;
use spenso::structure::representation::{IndexDisplay, IndexRow, RepName, RepresentationClass};
use spenso::structure::slot::IsAbstractSlot;
use symbolica::atom::{Atom, AtomView, FunctionBuilder};

const STYLE: &str = include_str!("../../typst/metadata.css");
static NEXT_INSPECTOR: AtomicU64 = AtomicU64::new(0);

fn card(kind: &str, body: &str) -> String {
    format!(
        "<style>{STYLE}</style><div data-spenso-metadata><div class=\"sm-type\"><code>{kind}</code></div>{body}</div>"
    )
}

pub(crate) fn representation_name(name: &SpensoRepresentationName) -> String {
    card(
        "RepresentationName",
        &format!(
            "<code class=\"sm-name\">{}</code><div class=\"sm-info\">{} · metric {}</div>",
            display::escape_html(name.rep.symbol().get_stripped_name()),
            name.duality(),
            name.metric_rule(),
        ),
    )
}

pub(super) fn metric(rep: Representation<LibraryRep>) -> String {
    let name = SpensoRepresentationName { rep: rep.rep };
    if !rep.rep.is_self_dual() {
        return "dual pairing".into();
    }
    let Dimension::Concrete(dimension) = rep.dim else {
        return name.metric_rule().into();
    };
    let mut signs = (0..dimension.min(16))
        .map(|i| if rep.is_neg(i) { "−" } else { "+" })
        .collect::<Vec<_>>();
    if dimension > 16 {
        signs.push("…");
    }
    format!("({})", signs.join(", "))
}

pub(crate) fn representation(rep: Representation<LibraryRep>) -> String {
    let name = SpensoRepresentationName { rep: rep.rep };
    card(
        "Representation",
        &format!(
            "<div class=\"sm-row\"><code class=\"sm-name\">{}</code><span>Dimension {}</span></div><div class=\"sm-info\">{} · {}</div>",
            display::escape_html(name.rep.symbol().get_stripped_name()),
            display::escape_html(&rep.dim.to_string()),
            name.duality(),
            metric(rep),
        ),
    )
}

impl IndexAliases {
    /// Script placement shared by concrete Typst and HTML interfaces. Dualizable
    /// representations reverse their base row; self-dual spaces retain it.
    pub(super) fn port_row(slot: PartialSlot) -> IndexRow {
        let representation = slot.rep_name();
        representation
            .metadata()
            .map(|metadata| {
                if metadata.class == RepresentationClass::Dualizable && representation.is_dual() {
                    metadata.index_row.opposite()
                } else {
                    metadata.index_row
                }
            })
            .unwrap_or(IndexRow::Top)
    }

    pub(super) fn port_display(&self, slot: PartialSlot) -> (Atom, Option<IndexDisplay>) {
        let port = composition::port_atom(slot);
        let mut matcher = SlotMatcher::default();
        let SlotMatch::Explicit(view) = matcher.classify(port.as_view()) else {
            return (port, None);
        };
        let index = view.index();
        let key = (slot.rep_name().symbol(), index.to_owned());
        if let Some((_, display)) = self.entries.get(&key) {
            return (index.to_owned(), Some(display.clone()));
        }
        if let Some(source) = self.source_indices.get(&key) {
            return (source.clone(), None);
        }
        let display = usize::try_from(index)
            .ok()
            .and_then(|position| slot.rep_name().metadata()?.index_palette.resolve(position))
            .or_else(|| match index {
                AtomView::Var(var) => IndexDisplay::from_symbol(var.get_symbol()),
                _ => None,
            });
        (index.to_owned(), display)
    }

    /// Shared with the component explorer's axis labels.
    pub(super) fn port_label(&self, slot: PartialSlot, settings: &DisplaySettings) -> String {
        let (index, display) = self.port_display(slot);
        display
            .map(|label| label.to_native_string())
            .unwrap_or_else(|| {
                display::format_atom_with_settings(&index, TensorDisplayMode::Plain, settings)
            })
    }

    pub(super) fn port_typst(&self, slot: PartialSlot, settings: &DisplaySettings) -> String {
        if matches!(slot.aind, PartialIndex::Open(_)) {
            return "?".into();
        }
        let (index, display) = self.port_display(slot);
        display
            .map(|display| display.to_typst_source())
            .unwrap_or_else(|| {
                display::format_atom_with_settings(&index, TensorDisplayMode::Typst, settings)
            })
    }

    pub(super) fn port_html(
        &self,
        py: Python<'_>,
        slot: PartialSlot,
        settings: &DisplaySettings,
    ) -> PyResult<String> {
        if matches!(slot.aind, PartialIndex::Open(_)) {
            return Ok("?".into());
        }
        let (index, display) = self.port_display(slot);
        if let Some(body) = display.and_then(|display| display.to_mathml()) {
            return Ok(format!(
                "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">{body}</math>"
            ));
        }
        paging::atom_to_html(py, &index, settings)
    }
}

pub(crate) fn slot(
    py: Python<'_>,
    slot: Slot<LibraryRep>,
    settings: &DisplaySettings,
) -> PyResult<String> {
    let atom = slot.to_atom();
    let aliases = IndexAliases::for_atom(&atom, &settings.index_style);
    let label = aliases.port_html(
        py,
        slot.rep().slot(PartialIndex::Explicit(slot.aind)),
        settings,
    )?;
    Ok(card(
        "Slot",
        &format!(
            "<div class=\"sm-row\"><span class=\"sm-index\">{}</span><code>{}</code></div><div class=\"sm-info\">Index <code>{}</code></div>",
            label,
            display::escape_html(&slot.rep().to_symbolic([]).to_string()),
            display::escape_html(&Atom::from(slot.aind).to_string()),
        ),
    ))
}

pub(crate) fn structure(
    py: Python<'_>,
    structure: &SpensoTensorStructure,
    settings: &DisplaySettings,
) -> PyResult<String> {
    let id = NEXT_INSPECTOR.fetch_add(1, Ordering::Relaxed);
    let selector = format!("[data-spenso-structure=\"{id}\"]");
    let slots = structure.interface.external_structure();
    let mut bundle = FunctionBuilder::new(AIND_SYMBOLS.aind);
    for slot in &slots {
        bundle = bundle.add_arg(composition::port_atom(*slot));
    }
    let aliases = IndexAliases::for_atom(&bundle.finish(), &settings.index_style);
    let mut tiles = String::new();
    let mut details = String::new();
    let mut rules = String::new();
    let mut count = 0;
    let mut add_tile = |caption: &str, preview: String, description: &str, detail: String| {
        rules.push_str(&format!("{selector}:has(input[value=\"{count}\"]:checked) .sm-detail[data-choice=\"{count}\"]{{display:block}}"));
        tiles.push_str(&format!(
            "<label class=\"sm-tile\"><input type=\"radio\" name=\"spenso-structure-{id}\" value=\"{count}\" aria-label=\"{}\" {}><span>{}</span><span class=\"sm-index\">{preview}</span><code>{}</code></label>",
            display::escape_html(caption), if count == 0 { "checked" } else { "" }, display::escape_html(caption), display::escape_html(description),
        ));
        details.push_str(&format!(
            "<div class=\"sm-detail\" data-choice=\"{count}\">{detail}</div>"
        ));
        count += 1;
    };
    for (position, slot) in slots.iter().enumerate() {
        let (kind, index) = match slot.aind {
            PartialIndex::Explicit(index) => (
                "Slot",
                format!(
                    "Index <code>{}</code>",
                    display::escape_html(&Atom::from(index).to_string())
                ),
            ),
            PartialIndex::Open(_) => ("Representation", "Unresolved index".to_owned()),
        };
        add_tile(
            &format!("{kind} {position}"),
            aliases.port_html(py, *slot, settings)?,
            &slot.rep().to_symbolic([]).to_string(),
            format!("<span class=\"sm-field\">{kind} {position}</span>{index}"),
        );
    }
    if count == 0 {
        tiles.push_str("<span class=\"sm-info\">Scalar · no external slots</span>");
    }
    let shape = slots
        .iter()
        .map(|slot| slot.rep().dim.to_string())
        .collect::<Vec<_>>()
        .join(" × ");
    Ok(format!(
        "<style>{STYLE}{rules}</style><div data-spenso-metadata data-spenso-structure=\"{id}\"><div class=\"sm-type\"><code>TensorStructure</code><span>Rank {}{}</span></div><div class=\"sm-ports\">{tiles}</div>{details}</div>",
        structure.interface.order(),
        if shape.is_empty() {
            String::new()
        } else {
            format!(" · {}", display::escape_html(&shape))
        },
    ))
}
