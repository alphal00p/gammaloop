use std::{
    collections::HashMap,
    sync::{LazyLock, RwLock},
};

use symbolica::atom::{AtomView, NamespacedSymbol, Symbol, SymbolBuilder};

use crate::{EntityKind, Model, ModelError};

struct ParameterLabel {
    latex: String,
    typst: String,
}

// Symbols have process-wide identities. Read labels at print time so importing
// a model again updates its presentation without redefining its symbols.
static LABELS: LazyLock<RwLock<HashMap<Symbol, ParameterLabel>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));

impl Model {
    /// Declare a parameter's printer before parsing model expressions.
    ///
    /// Import adapters which parse expressions before constructing a `Model`
    /// must call this for every parameter first, including forward references.
    /// Repeated declarations update labels; missing/empty labels restore the
    /// ordinary spelling. Existing user-declared symbol printers are preserved.
    pub fn register_parameter_symbol(name: &str, texname: Option<&str>) -> Result<(), ModelError> {
        let symbol_name = NamespacedSymbol::parse(&format!("UFO::{name}"));
        let symbol = if let Some(symbol) = Symbol::get_symbol(symbol_name.clone()) {
            // Existing user declarations (including their printers) are immutable.
            symbol
        } else {
            SymbolBuilder::new(symbol_name)
                .with_print_function(|view, options, _| {
                    if !options.mode.is_latex() && !options.mode.is_typst() {
                        return None;
                    }
                    let AtomView::Var(variable) = view else {
                        return None;
                    };
                    let labels = LABELS.read().unwrap();
                    let label = labels.get(&variable.get_symbol())?;
                    Some(if options.mode.is_latex() {
                        label.latex.clone()
                    } else {
                        label.typst.clone()
                    })
                })
                .build()
                .map_err(|message| ModelError::SymbolicParse {
                    kind: EntityKind::Parameter,
                    name: name.to_owned(),
                    field: "name",
                    expression: name.to_owned(),
                    message: message.to_string(),
                })?
        };
        let mut labels = LABELS.write().unwrap();
        if let Some(latex) = texname.filter(|name| !name.trim().is_empty()) {
            labels.insert(
                symbol,
                ParameterLabel {
                    latex: latex.to_owned(),
                    // MiTeX also handles labels with commands and grouped scripts;
                    // inserting raw LaTeX into Typst math would misrender those.
                    typst: format!("#{{ import \"@preview/mitex:0.2.6\": mi; mi({latex:?}) }}"),
                },
            );
        } else {
            labels.remove(&symbol);
        }
        Ok(())
    }
}
