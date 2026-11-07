use symbolica_notation::IndexDisplay;

use spenso::structure::partial::PartialIndex;

use crate::display;
use crate::display::{DisplaySettings, IndexAliases};
use idenso::tensor::composition;
use pyo3::{PyResult, Python};
use spenso::structure::partial::PartialSlot;
use spenso::structure::representation::{IndexRow, RepName, RepresentationClass};
use spenso::structure::slot::IsAbstractSlot;
use spenso::structure::slot::{SlotMatch, SlotMatcher};
use symbolica::atom::{Atom, AtomView};

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
    #[cfg(test)]
    pub(super) fn port_label(&self, slot: PartialSlot, settings: &DisplaySettings) -> String {
        let (index, display) = self.port_display(slot);
        display
            .map(|label| label.to_native_string())
            .unwrap_or_else(|| {
                display::format_atom_with_settings(
                    &index,
                    display::TensorDisplayMode::Plain,
                    settings,
                )
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
        display::atom_to_html(py, &index, settings, None)
    }
}
