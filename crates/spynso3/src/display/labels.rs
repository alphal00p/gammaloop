use std::{
    collections::BTreeMap,
    sync::{LazyLock, RwLock},
};

use spenso::network::tags::TENSOR_PRINT_HEAD_PREFIX;
use symbolica::atom::{Atom, AtomCore, AtomView, NamespacedSymbol, Symbol, SymbolBuilder};

use super::DisplaySettings;

struct LatexLabel {
    presentation: Symbol,
    typst: String,
}

static LABELS: LazyLock<RwLock<BTreeMap<Symbol, LatexLabel>>> =
    LazyLock::new(|| RwLock::new(BTreeMap::new()));

impl DisplaySettings {
    /// Register presentation-only LaTeX and native Typst names for an existing scalar symbol.
    ///
    /// Model importers can supply their labels even when the algebraic symbols
    /// already exist. Reloading a label updates display only; plain text and
    /// exact Atom payloads retain the original symbol identities.
    /// An empty label removes a previous registration.
    pub fn register_names(symbol: Symbol, latex: &str, typst: &str) {
        if latex.trim().is_empty() && typst.trim().is_empty() {
            LABELS.write().unwrap().remove(&symbol);
            return;
        }
        let ordinary = symbol.get_name();
        let ordinary = ordinary.rsplit("::").next().unwrap();
        let latex = if latex.trim().is_empty() {
            ordinary
        } else {
            latex
        };
        let typst = if typst.trim().is_empty() {
            format!("{ordinary:?}")
        } else {
            typst.to_owned()
        };
        // Include the original symbol: distinct parameters with identical labels
        // must not combine while constructing a temporary presentation Atom.
        let identity = symbol
            .get_name()
            .bytes()
            .chain([0])
            .chain(latex.bytes())
            .chain([0])
            .chain(typst.bytes())
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>();
        let name = NamespacedSymbol::parse(&format!("spenso::latex_label_{identity}"));
        let presentation = Symbol::get_symbol(name.clone()).unwrap_or_else(|| {
            let latex = latex.to_owned();
            let typst = typst.clone();
            SymbolBuilder::new(name)
                .with_tags([
                    format!("{TENSOR_PRINT_HEAD_PREFIX}latex:{latex}"),
                    format!("{TENSOR_PRINT_HEAD_PREFIX}typst:{typst}"),
                ])
                .with_print_function(move |_, options, _| {
                    if options.mode.is_latex() {
                        Some(latex.clone())
                    } else if options.mode.is_typst() {
                        Some(typst.clone())
                    } else {
                        None
                    }
                })
                .build()
                .expect("presentation names encode their complete declaration")
        });
        LABELS.write().unwrap().insert(
            symbol,
            LatexLabel {
                presentation,
                typst,
            },
        );
    }
}

pub(super) fn presentation_atom(atom: &Atom) -> Atom {
    let labels = LABELS.read().unwrap();
    if labels.is_empty() {
        return atom.clone();
    }
    atom.replace_map_bottom_up(|view, _, output| {
        if let AtomView::Var(variable) = view
            && let Some(label) = labels.get(&variable.get_symbol())
        {
            **output = Atom::var(label.presentation);
        }
    })
}

pub(super) fn typst_heads(atom: &Atom) -> BTreeMap<String, String> {
    let labels = LABELS.read().unwrap();
    atom.get_all_symbols(true)
        .into_iter()
        .filter_map(|symbol| {
            labels
                .get(&symbol)
                .map(|label| (symbol.get_name().to_owned(), label.typst.clone()))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use symbolica::symbol;

    use super::*;
    use crate::display::{TensorDisplayMode, atom_to_latex, format_atom_with_mode};

    #[test]
    fn latex_labels_preserve_distinct_symbols_and_plain_output() {
        let first = symbol!("latex_labels::first");
        let second = symbol!("latex_labels::second");
        let atom = Atom::var(first) + Atom::var(second);
        let plain = format_atom_with_mode(&atom, TensorDisplayMode::Plain, false);
        for symbol in [first, second] {
            DisplaySettings::register_names(symbol, r"\alpha_s", "alpha_s");
            DisplaySettings::register_names(symbol, r"\alpha_s", "alpha_s");
        }
        let source = format_atom_with_mode(&atom, TensorDisplayMode::Typst, false);
        assert_eq!(source.matches("alpha_s").count(), 2);
        assert!(!source.contains("mitex"));
        assert_eq!(
            format_atom_with_mode(&atom, TensorDisplayMode::Plain, false),
            plain
        );
        assert_eq!(atom.get_all_symbols(true).len(), 2);
        assert_eq!(typst_heads(&atom).len(), 2);
        assert!(
            atom_to_latex(&atom, &Default::default(), Some(1)).starts_with(r"$$\begin{gathered}")
        );
        DisplaySettings::register_names(first, r"\beta", "beta");
        assert!(format_atom_with_mode(&atom, TensorDisplayMode::Latex, false).contains(r"\beta"));
        for symbol in [first, second] {
            DisplaySettings::register_names(symbol, "", "");
        }
        assert_eq!(presentation_atom(&atom), atom);
    }
}
